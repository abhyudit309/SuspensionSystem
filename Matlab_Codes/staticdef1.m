function F = staticdef1(x,k,m,h,w,g)
    L = sqrt(w^2+h^2);
    F = 2*k*(h-x)*(L/sqrt((h-x)^2+w^2)-1)-m*g;
end