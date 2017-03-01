% v = FullWindow(u, off1, on1, on2, off2, kind)
%
% Window function that is zero for u <= off1,
%  one for on1 <= u <= on2, and zero for u >= off2.
%
% Must have off1 <= on1 <= on2 <= off2, or unexpected results will come.
%
% If off1 = on1, or off2 = on2, there will be an abrupt cutoff (step
%  function-like) instead of a smooth cutoff there.
%
% u can be a matrix or a scalar.
%
function v = FullWindow(u, off1, on1, on2, off2, kind)
    if off1 == on1,
        v = u >= on1;
    else
        v = HalfWindow(u, off1, on1, kind);
    end
    if off2 == on2,
        v = v .* (u <= off2);
    else
        v = v .* HalfWindow(u, off2, on2, kind);
    end
end