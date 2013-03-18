function plfit_mars, p, x=x, y=y, err=err

     model = p[0]*x[0, *]^(p[1])*x[1, *]^(p[2])*x[2, *]^(p[3])
     ; does this hackines... work?
     marserr = sqrt(abs(y-model))
     return, (y-model)/marserr
end