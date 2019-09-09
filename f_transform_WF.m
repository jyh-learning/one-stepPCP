function W0=f_transform_WF(W,F)
W0=W;
W0(F>0)=1-(1-F(F>0)).*(1-W(F>0));
W0(F<0)=(1+F(F<0)).*(W(F<0));
