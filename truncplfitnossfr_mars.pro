function truncplfit_mars, p, x=x, y=y, err=err

     model = p[0]*x[0, *]^(p[1])*x[1, *]^(p[2])*exp((-1)*x[0, *]/(p[4]))
     marserr = sqrt(abs(y-model))
     return, (y-model)/marserr
end