function make_vector,min_v,max_v,n_bin,log=log

bin = {bound_min:0.d0,$
       bound_max:0.d0,$
       mean:0.d0,$
       mean_2d:0.d0 }

vector = REPLICATE( bin, n_bin )

if (keyword_set(log)) then log=1 else log=0


;; -- boundaries:
if (log eq 0) then begin
;   vector.bound_min = min_v+(findgen(n_bin+1))(0:n_bin-1)/float(n_bin)*(max_v-min_v)
   vector.bound_min = min_v+(findgen(n_bin))/float(n_bin)*(max_v-min_v)
   vector.bound_max = min_v+(findgen(n_bin)+1.)/float(n_bin)*(max_v-min_v)
endif

if (log eq 1) then begin
   vector.bound_min = min_v*10.^(  findgen(n_bin)   /float(n_bin)*alog10(max_v/min_v))
   vector.bound_max = min_v*10.^( (findgen(n_bin)+1)/float(n_bin)*alog10(max_v/min_v))
endif


;; -- mean 1D:
for i=0L,n_bin-1 do begin
   vector(i).mean = (vector(i).bound_min+vector(i).bound_max)/2.
endfor

;; -- mean 2D:
for i=0L,n_bin-1 do begin
   vector(i).mean_2d = 2./3.*(vector[i].bound_max^3-vector[i].bound_min^3)/(vector[i].bound_max^2-vector[i].bound_min^2)
endfor


return,vector
end
