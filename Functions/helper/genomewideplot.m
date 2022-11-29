function f = genomewideplot(chr, pos, score, threshold, highlight)
% INPUT:
% chr: chromosome number (1 entry per SNP)
% pos: chromosomal position (1 entry per SNP)
% score: -log10(P) (1 entry per SNP)
% threshold: P-value threshold (note: do not take the log here), a default is plotted at p = 5e-8
% highlight: index of SNPs to be highlighted with a red circle
      
  if nargin < 4
      threshold = 5e-8;
  end
  if nargin < 5
      highlight = [];
  end

  % Alternating shades for chromosomes.
  clr1 = [109 188 231]/255; %light blue
  clr2 = [21 65 118]/255; %dark blue
  clr3 = [1 0 0]; % for highlighting
  
  % Get the set of chromosomes represented by the SNPs.
  chrs = unique(chr);
  chrs = chrs(:)';
  
  % This is the current position in the plot.
  p0 = 0;

  % Repeat for each chromosome.
  if nargout > 0, f = figure; end
  hold on
  M = zeros(max(22,length(chrs)),1); % save mean position per chromosome (for xlabel ticks)
  for c = chrs
    
    % Plot the SNPs on the chromosome.
    is = find(chr == c);
    maxpos = max(pos(is));
    p0 = p0 - min(pos(is)) + 10000; % min(pos) is used to avoid a gap in front of chromosomes
    if rem(c,2) == 1 % alternating color per chromosome
      clr = clr1;
    else
      clr = clr2;
    end
    hl = is(ismember(is,highlight));
    scatter(p0 + pos(is),score(is), 20 ,clr, 'filled'); % plotting SNPs
    scatter(p0 + pos(hl),score(hl), 100 ,clr3, 'd'); % highlighting SNPs

    % Keep track of chromosome means for labling
    M(c) = p0 + maxpos/2;
     
    % Move to the next chromosome.
    p0 = p0 + maxpos;
  end

  % Add threshold
  if threshold ~= 0, yline(-log10(threshold)); end
  hold off
  
  % Format X axis
  xticks(nonzeros(M));
  xticklabels(num2cell(chrs));
  xlim([0,p0])

  % Format Y axis
  ylabel('-log_1_0(p)')
  ylimit = max(max(score), 12);
  ylim([0 ylimit])

  % Format figure
  x0=10;
  y0=10;
  width=2000;
  height=400;
  set(gcf,'position',[x0,y0,width,height])
  set(gca, 'FontSize', 12);
  
  
  
  