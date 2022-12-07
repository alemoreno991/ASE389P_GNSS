function [PRN] = generatePRN(i)
% Generate the original LFSR sequences
n = 10;
ciVec_f1 = [3 10]';
ciVec_f2 = [2 3 6 8 9 10]';
a0Vec = [1 1 1 1 1 1 1 1 1 1]';
lfsrSeq_f1 = generateLfsrSequence(n, ciVec_f1, a0Vec);
lfsrSeq_f2 = generateLfsrSequence(n, ciVec_f2, a0Vec);

% Shift and get PRN
code_delay_chips = [5 6 7 8 17 18 139 140 141 251 252 254 255 256 257 ... 
    258 469 470 471 472 473 474 509 512 513 514 515 516 859 860 861 862 ...
    863 950 947 948 950 ];

lfsrSeq_f2_shifted = [
    lfsrSeq_f2(end-code_delay_chips(i)+1:end);... 
    lfsrSeq_f2(1:end-code_delay_chips(i))...
];
PRN = mod( lfsrSeq_f1 + lfsrSeq_f2_shifted, 2);

% Check if the PRN complies to IS-GPS-200L
PRN_checkVtr = [ 1440 1620 1710 1744 1133 1455 1131 1454 1626 1504 1642 ...
    1750 1764 1772 1775 1776 1156 1467 1633 1715 1746 1763 1063 1706 ...
    1743 1761 1770 1774 1127 1453 1625 1712 1745 1713 1134 1456 1713 ];

if ~isequal(oct2poly(PRN_checkVtr(i))', PRN(1:10))
    error_msg = ['Error while generating PRN_' num2str(i,'%02d')];
    display(error_msg)
end

PRN = 2*PRN - 1;
end

