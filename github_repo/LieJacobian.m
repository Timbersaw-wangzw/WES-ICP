function J = LieJacobian(A,B,order)
%LIEDERIVHESSIAN Lie algebra derivation vector and Hessian Matrix
% FORMATION: log_vex(Aexp(x)B)
% order: first-order or second-order expansion
switch order
    case 'firse-order'
        if isrot(A)
            J=A;
        else
            J=Ad(A);
        end
    case 'second-order'
        if isrot(A)
            J=(eye(3)-0.5*skew(d))*A;
        else
            ad_vec=ad(A*B);
            J=(eye(6)-0.5*ad_vec+ad_vec*ad_vec/12)*Ad(A);
%             J=(eye(6)-0.5*ad_vec)*Ad(A);
        end
end

