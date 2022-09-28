function F = staticdef2(x,k,m,h,w,g,xs)
    L = sqrt(w^2+h^2);
    F = 2*k*(h-x)*(L/sqrt((h-x)^2+w^2)-1)+k*(x-xs)-m*g;
end