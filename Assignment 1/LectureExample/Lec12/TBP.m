function ff = TBP (t, xx, mu)

rr = xx(1:3) ; r = norm(rr) ;

vv = xx(4:6) ;

ff(1:3,1) = vv  ;

ff(4:6,1) = -mu/r^3 * rr ;

end