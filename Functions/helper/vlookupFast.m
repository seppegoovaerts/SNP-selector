function [ind12,ind21,raw] = vlookupFast(Array1,Array2,rows)
         if nargin<3, rows = false; end
         if rows
            [~,ind12] = ismember(Array1,Array2,'rows');
         else
            [~,ind12] = ismember(Array1,Array2);
         end
         raw = ind12;
         ind21 = find(ind12);
         ind12 = ind12(ind21);
end