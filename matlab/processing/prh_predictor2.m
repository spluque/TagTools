## Working Octave port
## Valid: Matlab, Octave
## markjohnson@st-andrews.ac.uk
## Last modified: 15 Nov 2019

## -*- texinfo -*-
## @deftypefn {Function File} @var{prh} = prh_predictor2(@var{P}, @var{A}, @var{fs})
## @deftypefnx {Function File} @var{prh} = prh_predictor2(@var{P}, @var{A}, @var{fs}, @var{MAXD})
## @deftypefnx {Function File} @var{prh} = prh_predictor2(@var{P}, @var{A})
## @deftypefnx {Function File} @var{prh} = prh_predictor2(@var{P}, @var{A}, @var{MAXD})
## Predict tag position on diving animal parameterized by p0, r0, and h0
##
## p0, r0, and h0 are the cannonical angles between the principal axes of
## the tag and the animal. The tag orientation on the animal can change with
## time and this function provides a way to estimate the orientation at the
## start and end of each suitable dive. The function critically assumes that
## the animal makes a sequence of short dives between respirations and that
## the animal remains upright (i.e., does not roll) during these shallow
## dives. See prh_predictor1 for a method more suitable to animals that rest
## horizontally at the surface. The function provides a graphical interface
## showing the estimated tag-to-animal orientation throughout the
## deployment. Follow the directions above the top panel of the figure to
## edit or delete an orientation estimate.
##
## @var{P} is a dive depth vector or sensor structure with units of m
## H2O. A is an acceleration matrix or sensor structure with columns [ax ay
## az]. Acceleration can be in any consistent unit, e.g., g or m/s^2, and
## must have the same sampling rate as P.
##
## @var{fs} is the sampling rate of the sensor data in Hz (samples per
## second). This is only needed if A and M are not sensor structures.
##
## @var{MAXD} is the optional maximum depth of near-surface dives. The
## default value is 10 m. This is used to find contiguous surface intervals
## suitable for analysis.
##
## Returns:
##
## PRH = [cue,p0,r0,h0,q] with a row for each dive edge analysed.
## cue is the time in second-since-tagon of the dive edge
## analysed. [p0,r0,h0] are the deduced tag orientation angles in
## radians. q is the quality indicator with a low value (e.g., <0.05)
## indicating that the data fit more consistently with the assumptions of
## the method.
## @end deftypefn

function PRH = prh_predictor2(P, A, fs, MAXD)

  MINSEG = 30;			# minimum surface segment length in seconds
  MAXSEG = 300;			# maximum surface segment length in seconds
  GAP = 5;			# keep at least 5s away from a dive edge
  PRH = [];

  if nargin < 3,
    help prh_predictor2
    return
  endif

  if isstruct(P),
    [P, A, fs] = sens2var(P, A, 'regular');
    if isempty(P), return, endif
  else
    if nargin < 3,
      help prh_predictor2
      return
    endif
  endif

  if fs >= 7.5,
    df = round(fs / 5);
    P = decdc(P, df);
    A = decdc(A,df) ;
    fs = fs/df ;
  endif

  A = A .* repmat(norm2(A) .^ (-1), 1, 3); # normalise A to 1 g
  v = depth_rate(P, fs, 0.2);		   # vertical speed

  if isempty(MAXD),
    MAXD = 10 ;
  endif

  MAXD = max(MAXD, 2);
  T = find_dives(P, fs, MAXD); # find dives more than MAXD from surface
  if isempty(T),
    msg = 'No dives deeper than %d found in dive profile - change MAXD\n';
    fprintf(msg, MAXD);
    return
  endif

  ## check for segments with sufficient variation in orientation
  V = zeros(size(S, 1), 1);
  for k=1:size(S, 1),
    ks = round(S(k, 1) * fs) + 1:round(S(k, 2) * fs);
    V(k) = norm(std(A(ks, :)));
  end
  thr = median(V) + 1.5 * iqr(V) * [-1 1];
  S = S(V > thr(1) & V < thr(2), :);

  # check if there is a segment before first dive and after last dive
  s1 = [max(T.start(1) - MAXSEG, 0), T.start(1)];
  se = [T.end(end), min(T.end(end) + MAXSEG, (length(P) - 1) / fs)];
  k = find(P(round(fs * s1(1)) + 1:round(fs * s1(2))) > MAXD, 1, 'last');
  if ~isempty(k),
    s1(1) = s1(1) + k / fs;
  endif
  k = find(P(round(fs * se(1)) + 1:round(fs * se(2))) > MAXD, 1);
  if ~isempty(k),
    se(2) = se(1) + (k - 1) / fs;
  endif
  S = [s1; [T.end(1:end - 1) T.start(2:end)]; se];
  S = S(find(diff(S, [], 2) > MINSEG), :);

  # break up long surfacing intervals
  while 1,
    k = find(diff(S, [], 2) > MAXSEG, 1);
    if isempty(k), break, endif
    S = [S(1:k - 1, :); ...
	 S(k, 1) + [0 MAXSEG]; ...
	 S(k, 1) + MAXSEG S(k, 2); ...
	 S(k + 1:end, :)];
  endwhile

  # check for segments with sufficient variation in orientation
  V = zeros(size(S, 1), 1);
  for k=1:size(S, 1),
    ks = round(S(k, 1) * fs) + 1:round(S(k, 2) * fs);
    V(k) = norm(std(A(ks, :)));
  endfor
  thr = median(V) + 1.5 * iqr(V) * [-1 1];
  S = S(V > thr(1) & V < thr(2), :);

  PRH = NaN(size(S, 1), 5);

  for k=1:size(S, 1),		# apply prh inference method on segments
   prh = applymethod2(A, v, fs, S(k, :));
   if isempty(prh), continue, endif
   PRH(k, :) = [mean(S(k, 1:2)) prh];
  endfor

  figure(1), clf
  plot_fig1(P, fs, PRH);

  while 1,                # user input to adjust results
    [gx, gy, butt] = ginput(1);
    [~, k] = min(abs(PRH(:, 1) - gx));
    if butt == 1,
      fprintf(' %d: %3.1f degrees\n', round(gx), gy);
    else
      switch butt
        case {'q', 'Q'}
          break;
	case 'e'
	  ss = plot_fig2(A, v, fs, S(k, :), PRH(k, :));
	  if isempty(ss),
	    z = [1:k - 1 k + 1:size(S, 1)];
            S = S(z, :);
            PRH = PRH(z, :);
	  else
	    S(k, :) = ss;
	    prh = applymethod2(A, v, fs, ss);
	    PRH(k, :) = [mean(ss) prh];
	  endif
 	  plot_fig1(P,fs,PRH) ;
        case 'x'
	  z = [1:k - 1 k + 1:size(S, 1)];
	  S = S(z, :);
          PRH = PRH(z, :);
	  plot_fig1(P, fs, PRH);
	case 'z'
	  xl = get(gca, 'XLim');
	  xl = gx + diff(xl) / 4 * [-1 1];
	  xl(1) = max(xl(1), 0);
	  xl(2) = min(xl(2), length(P) / fs);
	  subplot(311), set(gca, 'XLim', xl);
	  subplot(312), set(gca, 'XLim', xl);
	  subplot(313), set(gca, 'XLim', xl);
	case 'Z'
	  xl = get(gca, 'XLim');
	  xl = gx + diff(xl) * [-1 1];
	  xl(1) = max(xl(1), 0);
	  xl(2) = min(xl(2), length(P) / fs);
	  subplot(311), set(gca,'XLim', xl);
	  subplot(312), set(gca, 'XLim', xl);
	  subplot(313), set(gca, 'XLim', xl);
      endswitch
    endif
  endwhile
  return
endfunction


## For animals that do sequences of roll-free shallow dives. Chooses r0 and
## h0 to minimize the mean-squared y-axis acceleration in segment As and
## then chooses p0 for a mean pitch angle of 0.
##
## Break into eigen-axes: assuming that the motion is mostly planar, the
## eigenvalues of QQ will indicate how planar: the largest two eigenvalues
## describe the energy in the plane of motion; the smallest eigenvalue
## describes the energy in the invariant direction i.e., the axis of
## rotation
function prh = applymethod2(A, v, fs, ss)
  ks = round(ss(1) * fs) + 1:round(ss(2) * fs);
  As = A(ks,:);
  vs = v(ks);

  # energy ratio between plane-of-motion and axis of rotation
  QQ = As' * As;			# form outer product of acceleration
  if any(isnan(QQ))
    prh = [];
    return
  endif

  [V, D] = svd(QQ);
  ## If the inverse condition pow>~0.05, the motion in As is more
  ## three-dimensional than two-dimensional
  pow = D(3, 3) / D(2, 2);

  ## Axis of rotation to restore V to tag Y axis
  aa = acos([0 1 0] * V(:, 3));
  Phi = cross([0; 1; 0], V(:, 3)) / sin(aa);
  ## Skew matrix
  S = [0, -Phi(3), Phi(2); Phi(3), 0, -Phi(1); -Phi(2), Phi(1), 0];
  ## Generate rotation matrix for rotation of aa degrees around axis Phi
  Q = eye(3) + (1 - cos(aa)) * S * S - sin(aa) * S;
  am = mean(As) * Q';
  p0 = atan2(am(1), am(3));
  Q = euler2rotmat([p0 0 0]) * Q;
  prh = [asin(Q(3, 1)) atan2(Q(3, 2), Q(3, 3)) atan2(Q(2, 1), Q(1, 1))];

  aa = As * Q(2, :)';
  prh(4) = mean([pow, std(aa)]);

  ## Check that h0 is not 180 degrees out by checking that the regression
  ## between Aa(:,1) and depth_rate is negative.
  Q = euler2rotmat(prh(1:3));	# make final transformation matrix
  Aa = As * Q';			# animal frame acceleration for the segment
  pp = polyfit(Aa(:, 1), vs, 1);
  if pp(1) > 0,
    ## if incorrect, add/subtract 180 degrees
    prh(3) = rem(prh(3) - pi, 2 * pi);
  endif

  for k=2:3,	# by convention, constrain r0 and h0 to the interval -pi:pi
    if abs(prh(k)) > pi,
      prh(k) = prh(k) - sign(prh(k)) * 2 * pi;
    endif
  endfor
  return
endfunction


function plot_fig1(P, fs, PRH)
  figure(1)
  if isempty(get(gcf, 'Children')),
    xl = [0 length(P) / fs];
  else
    xl = get(gca, 'XLim');
  endif

  subplot(311)
  plot((1:length(P)) / fs, P), set(gca, 'YDir', 'reverse'), grid
  ylabel('Depth, m')
  set(gca,'XLim', xl, 'XTickLabel', []);
  title('type e to edit, x to delete, z or Z to zoom in/out, or q to quit')
  subplot(312)
  plot(PRH(:, 1), PRH(:, 2:4) * 180 / pi, '*-'), grid
  set(gca, 'XLim', xl, 'XTickLabel', [])
  ylabel('PRH, degrees')
  subplot(313)
  plot(PRH(:, 1), min(PRH(:, 5), 0.15), '*-'), grid
  set(gca, 'XLim', xl, 'YLim', [0 0.15])
  xlabel('time cue')
  ylabel('Quality')
  return
endfunction


function seg = plot_fig2(A, v, fs, seg, prh)
  YEXT = 12.5;           % vertical extent of accelerometry plots +/-m/s^2
  rctx = [0 1 1 0; 1 1 0 0];
  rcty = [0 0 1 1; 0 1 1 0];

  while 1,
    figure(2), clf
    xl = seg + [-30 30];
    Aw = A * euler2rotmat(prh(2:4))';
    subplot(211)
    plot((1:size(A, 1)) / fs, A * 9.81), set(gca, 'XLim', xl), grid
    ylabel('Tag frame A, m/s^2')
    set(gca, 'YLim', YEXT * [-1 1], 'XLim', xl, 'XTickLabel', []);
    hold on
    plot(diff(seg(1:2)) * rctx + seg(1), 0.9 * YEXT * (2 * rcty - 1), 'k');
    title('type 1 or 2 to change segments, x to erase, or q to quit')
    subplot(212)
    plot((1:size(A, 1)) / fs, Aw * 9.81), set(gca, 'XLim', xl), grid
    ylabel('Animal frame A, m/s^2')
    set(gca, 'YLim', YEXT * [-1 1], 'XLim', xl);
    hold on
    plot(diff(seg(1:2)) * rctx + seg(1), 0.9 * YEXT * (2 * rcty - 1), 'k');
    mess = sprintf('%s  p0=%3.1f  r0=%3.1f  h0=%3.1f  quality=%4.3f', ...
		   prh(1), prh(2:4) * 180 / pi, prh(5));
    title(mess, 'FontSize', 12);
    [gx, gy, butt] = ginput(1);
    gx = max(min(gx, size(A, 1) / fs), 0);
    if butt == '1',
      seg(1) = gx;
    elseif butt == '2',
      seg(2) = gx;
    else
      if butt == 'x',
        seg = [];
      endif
      figure(1)
      return
    endif
    prh(2:end) = applymethod2(A, v, fs, seg);
  endwhile
  return
endfunction
