function fr=FRsimpAdEx(values,I,w0,V0,bmin)
% computes the firing rate of the simpAdEx for a given constant current I

% INPUT:
% * values:     vector containg the model parameters
% * I:          constant current step in pA
% * w0:         initial w-value; if w0 is empty, the steady-state  value 
%               (wr) is used
% * V0:         initial V-value; if V0 is empty, the steady-state value 
%               (Vr) is used
% * bmin:       lower bound for the model parameter b 

% OUTPUT:
% * fr:         firing rate in Hz


    % define parameters and values
    warning off;
    [Cm,gL,EL,sf,Vup,tcw,~,b,Vr,Vth]=names(values);

    tau=Cm/gL;
    f=tau/tcw;
    X_Vth=f*(I+gL*sf-gL*(Vth-EL));                  
    w_end=-gL*(Vth-EL)+gL*sf+I-X_Vth;
    if isempty(w0)
        if b~=0
            w_r=w_end+b;
        else
            w_r=0;
        end
    else
        w_r=w0;
    end
    if (isempty(bmin) || bmin<0)
        bmin=0;
    end
    if isempty(V0)
        V_r=Vr;
    else
        V_r=V0;
    end
   
    % compute exact firing rate
    if (X_Vth<=0 || f>=1 || b<bmin || Cm<=0 || gL<=0 || tcw<=0 || sf<=0 )
        % constraints not satisfied or I less than the rheobase
        fr=0;
    else
        X_Vr=f*(I+gL*sf*exp((V_r-Vth)./sf)-gL*(V_r-EL));  
        wV_Vr=-gL*(V_r-EL)+gL*sf*exp((V_r-Vth)./sf)+I;
        w1=wV_Vr-X_Vr;
        w2=wV_Vr+X_Vr;

        % first and second regime
        t1=0; t2=0;
        if V_r >= Vth
            if w_r>=wV_Vr
                w_ref=-gL*(Vth-EL)+gL*sf+I+X_Vth;
                Vfix=IntPoints_simpAdEx(values,I,w_r,w_ref,1);
                i=0; j=w_r; k=0;
                if (isempty(Vfix) || length(Vfix)==1)
                    t1=integral(@(x) ISI_int_piecewise(x,I,values,i,j,k),V_r,Vth);
                    t2=0;
                else
                    Vlb=min(Vfix);
                    t1=integral(@(x) ISI_int_piecewise(x,I,values,i,j,k),V_r,Vlb);
                    m=f*gL-gL;
                    i=m; j=(1-f)*I-m*EL; k=(1-f)*gL*sf;
                    t2=integral(@(x) ISI_int_piecewise(x,I,values,i,j,k),Vlb,Vth);
                end
                w_stop=w_end;
                Vlb=Vth;
            else
                Vlb=V_r;
                w_stop=w_r;
            end            
        else
            if(w_r < w2 && w_r > w1)
                m=f*gL-gL;
                i=m; j=(1-f)*I-m*EL; k=(1-f)*gL*sf;
                t2=integral(@(x) ISI_int_piecewise(x,I,values,i,j,k),V_r,Vth);
                w_stop=w_end;
            else
                i=0; j=w_r; k=0;
                if (w_r <= w1)
                    ns=-1;
                    w_ref=-gL*(Vth-EL)+gL*sf+I-X_Vth;
                else
                    ns=1;
                    w_ref=-gL*(Vth-EL)+gL*sf+I+X_Vth; % w_end
                end
                Vfix=IntPoints_simpAdEx(values,I,w_r,w_ref,ns);
                if isempty(Vfix) || length(Vfix)==1
                    t1=integral(@(x) ISI_int_piecewise(x,I,values,i,j,k),V_r,Vth);
                    w_stop=w_r;
                else
                    Vlb=min(Vfix);
                    t1=integral(@(x) ISI_int_piecewise(x,I,values,i,j,k),V_r,Vlb);
                    m=f*gL-gL;
                    i=m; j=(1-f)*I-m*EL; k=(1-f)*gL*sf;
                    t2=integral(@(x) ISI_int_piecewise(x,I,values,i,j,k),Vlb,Vth);
                    w_stop=w_end;
                end
            end
            Vlb=Vth;
        end
        
        % third regime
        if Vlb>=Vup
            t3=0;
        else
            i=0; j=w_stop; k=0;
            t3=integral(@(x) ISI_int_piecewise(x,I,values,i,j,k),Vlb,Vup);
        end
        % exact interspike interval and fr:
        ISI=t1+t2+t3;
        fr=1000/ISI;
    end
end


% ########### functions #############

function f=ISI_int_piecewise(V,I,values,i,j,k)
    [Cm,gL,EL,sf,~,~,~,~,~,Vth]=names(values);

    F=(1./Cm).*(I-(i.*V+j+k.*exp((V-Vth)./sf))+gL.*sf.*exp((V-Vth)./sf)-gL.*(V-EL));
    f=1./F;
end

function Vfix=IntPoints_simpAdEx(values,I,w0,w_ref,i)

    [Cm,gL,EL,sf,Vup,tcw,~,~,~,Vth]=names(values);

    if w0>w_ref
        f=Cm./(gL*tcw);
        G=@(V) ((1+i*f)*(I-gL*(V-EL)+gL*sf*exp((V-Vth)./sf)) - w0);
        options = optimset('Display','off');

        % first fixed point
        lb=EL+(I-(w0/(1+i*f)))/gL-0.1;
        ub=EL+sf+(I-(w0/(1+i*f)))/gL;
        while(sign(G(ub))/sign(G(lb))==1)
            lb=EL+(I-(w0/(1+i*f)))/gL-numcor;
            ub=EL+sf+(I-(w0/(1+i*f)))/gL;
            numcor=numcor+1;
            if(numcor>1000)
                disp(['[Cm gL EL sf Vup tcw a b Vr Vth]=[' num2str([Cm,gL,EL,sf,Vup,tcw,a,b,Vr,Vth]) ']']);
                disp(['w0=' num2str(w0)]);
                disp(['I=' num2str(I)]);
                error('Numerical problems! Please check the provided input parameters and/or change the maximal iteration steps "numcor" (line 143)!');
            end
        end
        Vfix(1)=fzero(G,[lb ub],options);
        
        % second fixed point
        lb=EL+sf+(I-(w0/(1+i*f)))/gL;
        ub=Vup;
        numcor=0;
        while(sign(G(ub))/sign(G(lb))==1)
            lb=EL+sf+(I-(w0/(1+i*f)))/gL-numcor;
            ub=Vup;
            numcor=numcor+1;
            if(numcor>1000)
                disp(['[Cm gL EL sf Vup tcw a b Vr Vth]=[' num2str([Cm,gL,EL,sf,Vup,tcw,a,b,Vr,Vth]) ']']);
                disp(['w0=' num2str(w0)]);
                disp(['I=' num2str(I)]);
                error('Numerical problems! Please check the provided input parameters and/or change the maximal iteration steps "numcor" (line 143)!');
            end
        end
        Vfix(2)=fzero(G,[lb ub],options);
        
    elseif w0==w_ref
        Vfix=Vth;        
    else
        Vfix=[];        
    end 
    
end


% (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim
