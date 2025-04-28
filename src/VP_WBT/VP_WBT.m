classdef VP_WBT < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        figure1                      matlab.ui.Figure
        WeightSettingEditField       matlab.ui.control.EditField
        WeightSettingEditFieldLabel  matlab.ui.control.Label
        WeightDropDown               matlab.ui.control.DropDown
        WeightDropDownLabel          matlab.ui.control.Label
        edit16                       matlab.ui.control.TextArea
        TTextAreaLabel               matlab.ui.control.Label
        popupmenu7                   matlab.ui.control.DropDown
        edit15                       matlab.ui.control.EditField
        text30                       matlab.ui.control.Label
        edit14                       matlab.ui.control.EditField
        edit13                       matlab.ui.control.EditField
        text29                       matlab.ui.control.Label
        text28                       matlab.ui.control.Label
        pushbutton6                  matlab.ui.control.Button
        popupmenu6                   matlab.ui.control.DropDown
        popupmenu5                   matlab.ui.control.DropDown
        text27                       matlab.ui.control.Label
        popupmenu4                   matlab.ui.control.DropDown
        text26                       matlab.ui.control.Label
        text25                       matlab.ui.control.Label
        text24                       matlab.ui.control.Label
        text23                       matlab.ui.control.Label
        edit11                       matlab.ui.control.EditField
        text22                       matlab.ui.control.Label
        edit10                       matlab.ui.control.EditField
        text20                       matlab.ui.control.Label
        edit9                        matlab.ui.control.EditField
        text19                       matlab.ui.control.Label
        edit8                        matlab.ui.control.EditField
        edit7                        matlab.ui.control.EditField
        text16                       matlab.ui.control.Label
        edit6                        matlab.ui.control.EditField
        text15                       matlab.ui.control.Label
        edit5                        matlab.ui.control.EditField
        text13                       matlab.ui.control.Label
        text12                       matlab.ui.control.Label
        text11                       matlab.ui.control.Label
        pushbutton4                  matlab.ui.control.Button
        pushbutton3                  matlab.ui.control.Button
        popupmenu1                   matlab.ui.control.DropDown
        text10                       matlab.ui.control.Label
        text9                        matlab.ui.control.Label
        text8                        matlab.ui.control.Label
        text7                        matlab.ui.control.Label
        pushbutton1                  matlab.ui.control.Button
        edit4                        matlab.ui.control.EditField
        text6                        matlab.ui.control.Label
        text5                        matlab.ui.control.Label
        edit3                        matlab.ui.control.EditField
        edit2                        matlab.ui.control.EditField
        text4                        matlab.ui.control.Label
        text3                        matlab.ui.control.Label
        edit1                        matlab.ui.control.EditField
        axes1                        matlab.ui.control.UIAxes
    end

    
    
    methods (Access = private)
        function y=truncation(app,x,f)
            y=f(x);
            y(isinf(x))=0;
        end
        function a = Chebyshev(app, x, n)
            %代入n，计算Tn
            xsize=size(x,2);
            a=zeros(xsize,n+1);
            a(:,1)=1;
            a(:,2)=x.';
            for i=3:n+1
                a(:,i)=2*x.'.*a(:,i-1)-a(:,i-2);
            end
        end
        
        function ys = exptexting2(app, u, n, scale, s, v, stepsize, acc)
            %SOG
            gg=@(x) eval(s);
            g=@(x)truncation(app,x,gg);
            scale=n*scale;
            A=zeros(1,2*n);
            if acc<1
                for i=0:2*n-1
                    t=i;
                    q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
                    A(i+1)=quadgk(q,0,pi,'AbsTol',acc,'MaxIntervalCount',1e10,'RelTol',0);
                end
            else
                w=acc;
                ww=ceil(sqrt(w));
                d=2;
                S=2*ww/d;
                N=100;
                [b1,b2]=grule(app, N);
                for i=0:2*n-1
                    t=i;
                    q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
                    for tt=1:S
                        A(i+1)=A(i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
                    end
                    A(i+1)=A(i+1)/ww;
                end
            end
            A(1)=A(1)/2;
            u=(v+stepsize):stepsize:u;
            y=2*exp(-(u.^2)/scale)-1;
            for i=1:size(u,2)
                T=Chebyshev(app, y(i),2*n-1);
                K(i)=0;
                for k=0:n
                    K(i)=K(i)+A(k+1)*T(k+1);
                end
                for k=1:n-1
                    K(i)=K(i)+(1-k/n)*A(k+n+1)*T(k+n+1);
                end
            end
            z=g(u);
            y=abs(z-K);
            ys=[u;K;y];
        end
        
        function ys = exptexting3(app, u, n, scale, s, v, stepsize, acc)
            gg=@(x) eval(s);
            g=@(x)truncation(app,x,gg);
            scale=n*scale;
            A=zeros(1,2*n);
            if acc<1
                for i=0:2*n-1
                    t=i;
                    q=@(x) (2/pi)*cos(t*x).*g((-scale*log((1+cos(x))/2)));
                    A(i+1)=quadgk(q,0,pi,'AbsTol',acc,'MaxIntervalCount',1e10,'RelTol',0);
                end
            else
                w=acc;
                ww=ceil(sqrt(w));
                d=2;
                S=2*ww/d;
                N=100;
                [b1,b2]=grule(app, N);
                for i=0:2*n-1
                    t=i;
                    q=@(x) (2/pi)*cos(t*x).*g((-scale*log((1+cos(x))/2)));
                    for tt=1:S
                        A(i+1)=A(i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
                    end
                    A(i+1)=A(i+1)/ww;
                end
            end
            A(1)=A(1)/2;
            u=(v+stepsize):stepsize:u;
            y=2*exp(-(u)/scale)-1;
            for i=1:size(u,2)
                T=Chebyshev(app, y(i),2*n-1);
                K(i)=0;
                for k=0:n
                    K(i)=K(i)+A(k+1)*T(k+1);
                end
                for k=1:n-1
                    K(i)=K(i)+(1-k/n)*A(k+n+1)*T(k+n+1);
                end
            end
            z=g(u);
            y=abs(z-K);
            ys=[u;K;y];
        end
        
        function ys = exptexting5(app, u, n, scale, s1, s2, sp, v, stepsize, acc)
            g1=@(x) eval(s1);
            g2=@(x) eval(s2);
            scale=n*scale;sps=acos(2*exp(-sp*sp/scale)-1);
            A=zeros(1,2*n);
            if acc<1
                for i=0:2*n-1
                    t=i;
                    q1=@(x) (2/pi)*cos(t*x).*g1(sqrt(-scale*log((1+cos(x))/2)));
                    q2=@(x) (2/pi)*cos(t*x).*g2(sqrt(-scale*log((1+cos(x))/2)));
                    A(i+1)=quadgk(q1,0,sps,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0)+quadgk(q2,sps,pi,'AbsTol',1.0e-12,'MaxIntervalCount',1e9,'RelTol',0);
                end
            else
                w=acc;
                ww=ceil(sqrt(w));
                d=2;
                S=2*ww/d;
                N=100;
                [b1,b2]=grule(app, N);
                for i=0:2*n-1
                    t=i;
                    set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                    q1=@(x) (2/pi)*cos(t*x).*g1(sqrt(-scale*log((1+cos(x))/2)));
                    q2=@(x) (2/pi)*cos(t*x).*g2(sqrt(-scale*log((1+cos(x))/2)));
                    for tt=1:S
                        A(i+1)=A(i+1)+q1(((2*tt-1+b1-ww)/ww)*sps/2+sps/2)*b2'*sps/2;
                        A(i+1)=A(i+1)+q2(((2*tt-1+b1-ww)/ww)*(pi-sps)/2+(sps+pi)/2)*b2'*(pi-sps)/2;
                    end
                    A(i+1)=A(i+1)/ww;
                end
            end
            A(1)=A(1)/2;
            u=(v+stepsize):stepsize:u;
            NNN=size(u,2);
            y=2*exp(-(u.^2)/scale)-1;
            for i=1:size(u,2)
                T=Chebyshev(app, y(i),2*n-1);
                K(i)=0;
                for k=0:n
                    K(i)=K(i)+A(k+1)*T(k+1);
                end
                for k=1:n-1
                    K(i)=K(i)+(1-k/n)*A(k+n+1)*T(k+n+1);
                end
            end
            eta1=[ones(1,(sp-v)/stepsize) zeros(1,NNN-(sp-v)/stepsize)];
            eta2=[zeros(1,(sp-v)/stepsize) ones(1,NNN-(sp-v)/stepsize)];
            z=g1(u).*eta1+g2(u).*eta2;
            y=abs(z-K);
            ys=[u;K;y];
        end
        
        function ys = exptexting6(app, u, n, scale, s1, s2, sp, v, stepsize, acc)
            g1=@(x) eval(s1);
            g2=@(x) eval(s2);
            scale=n*scale;sps=acos(2*exp(-sp/scale)-1);
            A=zeros(1,2*n);
            if acc<1
                for i=0:2*n-1
                    t=i;
                    q1=@(x) (2/pi)*cos(t*x).*g1((-scale*log((1+cos(x))/2)));
                    q2=@(x) (2/pi)*cos(t*x).*g2((-scale*log((1+cos(x))/2)));
                    A(i+1)=quadgk(q1,0,sps,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0)+quadgk(q2,sps,pi,'AbsTol',1.0e-12,'MaxIntervalCount',1e9,'RelTol',0);
                end
            else
                w=acc;
                ww=ceil(sqrt(w));
                d=2;
                S=2*ww/d;
                N=100;
                [b1,b2]=grule(app, N);
                for i=0:2*n-1
                    t=i;
                    set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                    q1=@(x) (2/pi)*cos(t*x).*g1((-scale*log((1+cos(x))/2)));
                    q2=@(x) (2/pi)*cos(t*x).*g2((-scale*log((1+cos(x))/2)));
                    for tt=1:S
                        A(i+1)=A(i+1)+q1(((2*tt-1+b1-ww)/ww)*sps/2+sps/2)*b2'*sps/2;
                        A(i+1)=A(i+1)+q2(((2*tt-1+b1-ww)/ww)*(pi-sps)/2+(sps+pi)/2)*b2'*(pi-sps)/2;
                    end
                    A(i+1)=A(i+1)/ww;
                end
            end
            A(1)=A(1)/2;
            u=(v+stepsize):stepsize:u;
            NNN=size(u,2);
            y=2*exp(-(u)/scale)-1;
            for i=1:size(u,2)
                T=Chebyshev(app, y(i),2*n-1);
                K(i)=0;
                for k=0:n
                    K(i)=K(i)+A(k+1)*T(k+1);
                end
                for k=1:n-1
                    K(i)=K(i)+(1-k/n)*A(k+n+1)*T(k+n+1);
                end
            end
            eta1=[ones(1,(sp-v)/stepsize) zeros(1,NNN-(sp-v)/stepsize)];
            eta2=[zeros(1,(sp-v)/stepsize) ones(1,NNN-(sp-v)/stepsize)];
            z=g1(u).*eta1+g2(u).*eta2;
            y=abs(z-K);
            ys=[u;K;y];
        end
        
        function [bp, wf] = grule(app, n)
            % [bp,wf]=grule(n)
            %  This function computes Gauss base points and weight factors
            %  using the algorithm given by Davis and Rabinowitz in 'Methods
            %  of Numerical Integration', page 365, Academic Press, 1975.
            bp=zeros(n,1); wf=bp; iter=2; m=fix((n+1)/2); e1=n*(n+1);
            mm=4*m-1; t=(pi/(4*n+2))*(3:4:mm); nn=(1-(1-1/n)/(8*n*n));
            xo=nn*cos(t);
            for j=1:iter
               pkm1=1; pk=xo;
               for k=2:n
                  t1=xo.*pk; pkp1=t1-pkm1-(t1-pkm1)/k+t1;
                  pkm1=pk; pk=pkp1;
               end
               den=1.-xo.*xo; d1=n*(pkm1-xo.*pk); dpn=d1./den;
               d2pn=(2.*xo.*dpn-e1.*pk)./den;
               d3pn=(4*xo.*d2pn+(2-e1).*dpn)./den;
               d4pn=(6*xo.*d3pn+(6-e1).*d2pn)./den;
               u=pk./dpn; v=d2pn./dpn;
               h=-u.*(1+(.5*u).*(v+u.*(v.*v-u.*d3pn./(3*dpn))));
               p=pk+h.*(dpn+(.5*h).*(d2pn+(h/3).*(d3pn+.25*h.*d4pn)));
               dp=dpn+h.*(d2pn+(.5*h).*(d3pn+h.*d4pn/3));
               h=h-p./dp; xo=xo+h;
            end
            bp=-xo-h;
            fx=d1-h.*e1.*(pk+(h/2).*(dpn+(h/3).*(...
                d2pn+(h/4).*(d3pn+(.2*h).*d4pn))));
            wf=2*(1-bp.^2)./(fx.*fx);
            if (m+m) > n, bp(m)=0; end
            if ~((m+m) == n), m=m-1; end
            jj=1:m; n1j=(n+1-jj); bp(n1j)=-bp(jj); wf(n1j)=wf(jj);
        end
        
        function err = mexptexting2(app, u, n, scales, s, v, stepsize, acc)
            %SOG
            gg=@(x) eval(s);
            g=@(x)truncation(app,x,gg);
            scales=n*scales;
            sc=size(scales,2);
            A=zeros(sc,2*n);
            if acc<1
                for j=1:sc
                    scale=scales(j);
                    for i=0:2*n-1
                        t=i;
                        q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
                        A(j,i+1)=quadgk(q,0,pi,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0);
                    end
                end
            else
                 w=acc;
                ww=ceil(sqrt(w));
                d=2;
                S=2*ww/d;
                N=100;
                [b1,b2]=grule(app, N);
                for j=1:sc
                    scale=scales(j);
                    for i=0:2*n-1
                        t=i;
                        q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
                        for tt=1:S
                            A(j,i+1)=A(j,i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
                        end
                        A(j,i+1)=A(j,i+1)/ww;
                    end
                end
            end
            A(:,1)=A(:,1)/2;
            u=(v+stepsize):stepsize:u;
            scales=scales.';
            y=2*exp(-(u.^2)./scales)-1;
            K=zeros(sc,size(u,2));
            for i=1:size(u,2)
                T=Chebyshev(app, y(:,i).',2*n-1);
                for k=0:n
                    K(:,i)=K(:,i)+A(:,k+1).*T(:,k+1);
                end
                for k=1:n-1
                    K(:,i)=K(:,i)+(1-k/n)*A(:,k+n+1).*T(:,k+n+1);
                end
            end
            z=g(u).*ones(sc,1);
            y=abs(z-K);
            err=max(y.');
        end
        
        function err = mexptexting3(app, u, n, scales, s, v, stepsize, acc)
            %SOG
            g=@(x) eval(s);
            scales=n*scales;
            sc=size(scales,2);
            A=zeros(sc,2*n);
            if acc<1
                for j=1:sc
                    scale=scales(j);
                    for i=0:2*n-1
                        t=i;
                        q=@(x) (2/pi)*cos(t*x).*g((-scale*log((1+cos(x))/2)));
                        A(j,i+1)=quadgk(q,0,pi,'AbsTol',acc,'MaxIntervalCount',1e9,'RelTol',0);
                    end
                end
            else
                 w=acc;
                ww=ceil(sqrt(w));
                d=2;
                S=2*ww/d;
                N=100;
                [b1,b2]=grule(app, N);
                for j=1:sc
                    scale=scales(j);
                    for i=0:2*n-1
                        t=i;
                        q=@(x) (2/pi)*cos(t*x).*g((-scale*log((1+cos(x))/2)));
                        for tt=1:S
                            A(j,i+1)=A(j,i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
                        end
                        A(j,i+1)=A(j,i+1)/ww;
                    end
                end
            end
            A(:,1)=A(:,1)/2;
            u=(v+stepsize):stepsize:u;
            scales=scales.';
            y=2*exp(-(u)./scales)-1;
            K=zeros(sc,size(u,2));
            for i=1:size(u,2)
                T=Chebyshev(app, y(:,i).',2*n-1);
                for k=0:n
                    K(:,i)=K(:,i)+A(:,k+1).*T(:,k+1);
                end
                for k=1:n-1
                    K(:,i)=K(:,i)+(1-k/n)*A(:,k+n+1).*T(:,k+n+1);
                end
            end
            z=g(u).*ones(sc,1);
            y=abs(z-K);
            err=max(y.');
        end
        
        function popupmenu2_CreateFcn(app, hObject, eventdata, handles)
            % --- Executes during object creation, after setting all properties.
            
            % hObject    handle to popupmenu2 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    empty - handles not created until after all CreateFcns called
            
            % Hint: popupmenu controls usually have a white background on Windows.
            %       See ISPC and COMPUTER.
            if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
                set(hObject,'BackgroundColor','white');
            end
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function VP_WBT_OpeningFcn(app)
                        % --- Executes just before VPMR is made visible.
            
            % Ensure that the app appears on screen when run
            movegui(app.figure1, 'onscreen');
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app); %#ok<ASGLU>
            
            % This function has no output args, see OutputFcn.
            % hObject    handle to figure
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            % varargin   command line arguments to VPMR (see VARARGIN)
            
            % Choose default command line output for VPMR
            handles.output = hObject;
            handles.popindex=1;
            handles.err=[];
            handles.XX=0;
            handles.pp=[];
            handles.ww=[];
            handles.Anew=[];
            handles.Bnew=[];
            handles.Cnew=[];
            handles.radioindex=0;
            handles.radioindex2=0;
            handles.digitss=300;
            handles.xmax=0;
            handles.xmin=0;
            handles.n=0;
            handles.stepsize=0;
            handles.sclae=0;
            handles.errors=0;
            handles.errs=[];
            handles.s=[];
            handles.weight =[];
            handles.weightset =[];
            handles.T = 0;
            handles.accindex=1;
            handles.mrterms=0;
            handles.intindex=2;
            handles.records={'Record'};
            set(handles.pushbutton1,'Visible','Off');
            set(handles.pushbutton3,'Visible','Off');
            set(handles.pushbutton4,'Visible','Off');
            set(handles.pushbutton6,'Visible','Off');
            set(handles.text3,'Visible','Off');
            set(handles.text4,'Visible','Off');
            set(handles.text5,'Visible','Off');
            set(handles.edit1,'Visible','Off');
            set(handles.edit2,'Visible','Off');
            set(handles.edit3,'Visible','Off');
            set(handles.text7,'Visible','Off');
            set(handles.text8,'Visible','Off');
            set(handles.text11,'Visible','Off');
            set(handles.text13,'Visible','Off');
            set(handles.edit5,'Visible','Off');
            set(handles.text15,'Visible','Off');
            set(handles.edit6,'Visible','Off');
            set(handles.text16,'Visible','Off');
            set(handles.text19,'Visible','Off');
            set(handles.text20,'Visible','Off');
            set(handles.edit7,'Visible','Off');
            set(handles.edit8,'Visible','Off');
            set(handles.edit9,'Visible','Off');
            set(handles.edit10,'Visible','Off');
            set(handles.edit11,'Visible','Off');
            set(handles.text22,'Visible','Off');
            set(handles.text23,'Visible','Off');
            set(handles.text24,'Visible','Off');
            set(handles.text25,'Visible','Off');
            set(handles.text26,'Visible','Off');
            set(handles.popupmenu4,'Visible','Off');
            set(handles.popupmenu5,'Visible','Off');
            set(handles.popupmenu6,'Visible','Off');
            set(handles.text27,'Visible','Off');
            set(handles.text27,'Visible','Off');
            set(handles.edit13,'Visible','Off');
            set(handles.edit14,'Visible','Off');
            set(handles.edit15,'Visible','Off');
            set(handles.text28,'Visible','Off');
            set(handles.text29,'Visible','Off');
            set(handles.text30,'Visible','Off');
            set(app.edit16,'Visible','On');
            set(app.WeightDropDown,'Visible','Off');
            set(app.WeightDropDownLabel,'Visible','Off')
            set(app.WeightSettingEditField,'Visible','Off')
            set(app.WeightSettingEditFieldLabel,'Visible','Off')
            % Update handles structure
            guidata(hObject, handles);
        end

        % Value changed function: popupmenu1
        function popupmenu1_Callback(app, event)
            % --- Executes on selection change in popupmenu1.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to popupmenu1 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles.popindex=get(handles.popupmenu1,'Value');
            if handles.popindex==2
                set(handles.pushbutton1,'Visible','On');
                set(handles.pushbutton3,'Visible','On');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','On');
                set(handles.text4,'Visible','On');
                set(handles.text5,'Visible','On');
                set(handles.edit1,'Visible','On');
                set(handles.edit2,'Visible','On');
                set(handles.edit3,'Visible','On');
                set(handles.edit5,'Visible','Off');
                set(app.edit16,'Visible','On');
                set(handles.text7,'Visible','On');
                set(handles.text8,'Visible','On');
                set(handles.text11,'Visible','Off');
                set(handles.text12,'Visible','Off');
                set(handles.text13,'Visible','Off');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','Off');
                set(handles.edit7,'Visible','Off');
                set(handles.edit8,'Visible','On');
                set(handles.edit9,'Visible','On');
                set(handles.text19,'Visible','On');
                set(handles.text20,'Visible','On');
                set(handles.text23,'Visible','Off');
                set(handles.text24,'Visible','Off');
                set(handles.text25,'Visible','Off');
                set(handles.text26,'Visible','On');
                set(handles.popupmenu4,'Visible','On');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','Off');
                set(handles.edit14,'Visible','Off');
                set(handles.edit15,'Visible','Off');
                set(handles.text28,'Visible','Off');
                set(handles.text29,'Visible','Off');
                set(handles.text30,'Visible','Off');
            elseif handles.popindex==5
                set(handles.pushbutton1,'Visible','On');
                set(handles.pushbutton3,'Visible','On');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','On');
                set(handles.text4,'Visible','On');
                set(handles.text5,'Visible','On');
                set(handles.edit1,'Visible','On');
                set(handles.edit2,'Visible','On');
                set(handles.edit3,'Visible','On');
                set(handles.edit5,'Visible','Off');
                set(handles.text7,'Visible','On');
                set(handles.text8,'Visible','On');
                set(handles.text11,'Visible','On');
                set(handles.text12,'Visible','On');
                set(handles.text13,'Visible','Off');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','On');
                set(handles.edit7,'Visible','On');
                set(handles.edit8,'Visible','On');
                set(handles.edit9,'Visible','On');
                set(handles.text19,'Visible','On');
                set(handles.text20,'Visible','On');
                set(handles.text23,'Visible','Off');
                set(handles.text24,'Visible','Off');
                set(handles.text25,'Visible','Off');
                set(handles.text26,'Visible','On');
                set(handles.popupmenu4,'Visible','On');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','Off');
                set(handles.edit14,'Visible','Off');
                set(handles.edit15,'Visible','Off');
                set(handles.text28,'Visible','Off');
                set(handles.text29,'Visible','Off');
                set(handles.text30,'Visible','Off');
                set(app.WeightDropDown,'Visible','On');
                set(app.WeightDropDownLabel,'Visible','On')
                set(app.WeightSettingEditField,'Visible','On')
                set(app.WeightSettingEditFieldLabel,'Visible','On')
            elseif handles.popindex==3
                set(handles.pushbutton1,'Visible','On');
                set(handles.pushbutton3,'Visible','On');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','On');
                set(handles.text4,'Visible','Off');
                set(handles.text5,'Visible','On');
                set(handles.edit1,'Visible','On');
                set(handles.edit2,'Visible','Off');
                set(handles.edit3,'Visible','On');
                set(handles.edit5,'Visible','Off');
                set(handles.text7,'Visible','On');
                set(handles.text8,'Visible','On');
                set(handles.text11,'Visible','On');
                set(handles.text12,'Visible','On');
                set(handles.text13,'Visible','Off');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','Off');
                set(handles.edit7,'Visible','Off');
                set(handles.edit8,'Visible','On');
                set(handles.edit9,'Visible','On');
                set(handles.text19,'Visible','On');
                set(handles.text20,'Visible','On');
                set(handles.text23,'Visible','Off');
                set(handles.text24,'Visible','Off');
                set(handles.text25,'Visible','Off');
                set(handles.text26,'Visible','On');
                set(handles.popupmenu4,'Visible','On');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','On');
                set(handles.edit14,'Visible','On');
                set(handles.edit15,'Visible','On');
                set(handles.text28,'Visible','On');
                set(handles.text29,'Visible','On');
                set(handles.text30,'Visible','On');
            elseif handles.popindex==4
                set(handles.pushbutton1,'Visible','On');
                set(handles.pushbutton3,'Visible','On');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','On');
                set(handles.text4,'Visible','On');
                set(handles.text5,'Visible','On');
                set(handles.text7,'Visible','On');
                set(handles.text8,'Visible','On');
                set(handles.edit1,'Visible','On');
                set(handles.edit2,'Visible','On');
                set(handles.edit3,'Visible','On');
                set(handles.edit5,'Visible','On');
                set(handles.text11,'Visible','On');
                set(handles.text12,'Visible','On');
                set(handles.text13,'Visible','On');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','On');
                set(handles.edit7,'Visible','On');
                set(handles.edit8,'Visible','On');
                set(handles.edit9,'Visible','On');
                set(handles.text19,'Visible','On');
                set(handles.text20,'Visible','On');
                set(handles.text23,'Visible','Off');
                set(handles.text24,'Visible','Off');
                set(handles.text25,'Visible','Off');
                set(handles.text26,'Visible','On');
                set(handles.popupmenu4,'Visible','On');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','Off');
                set(handles.edit14,'Visible','Off');
                set(handles.edit15,'Visible','Off');
                set(handles.text28,'Visible','Off');
                set(handles.text29,'Visible','Off');
                set(handles.text30,'Visible','Off');
                set(app.WeightDropDown,'Visible','On');
                set(app.WeightDropDownLabel,'Visible','On')
                set(app.WeightSettingEditField,'Visible','On')
                set(app.WeightSettingEditFieldLabel,'Visible','On')
            elseif handles.popindex==1
                set(handles.pushbutton1,'Visible','Off');
                set(handles.pushbutton3,'Visible','Off');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','Off');
                set(handles.text4,'Visible','Off');
                set(handles.text5,'Visible','Off');
                set(handles.edit1,'Visible','Off');
                set(handles.edit2,'Visible','Off');
                set(handles.edit3,'Visible','Off');
                set(handles.edit5,'Visible','Off');
                set(handles.text7,'Visible','Off');
                set(handles.text8,'Visible','Off');
                set(handles.text11,'Visible','Off');
                set(handles.text12,'Visible','Off');
                set(handles.text13,'Visible','Off');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','Off');
                set(handles.edit7,'Visible','Off');
                set(handles.edit8,'Visible','Off');
                set(handles.edit9,'Visible','Off');
                set(handles.text19,'Visible','Off');
                set(handles.text20,'Visible','Off');
                set(handles.text23,'Visible','Off');
                set(handles.text24,'Visible','Off');
                set(handles.text25,'Visible','Off');
                set(handles.text26,'Visible','Off');
                set(handles.popupmenu4,'Visible','Off');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','Off');
                set(handles.edit14,'Visible','Off');
                set(handles.edit15,'Visible','Off');
                set(handles.text28,'Visible','Off');
                set(handles.text29,'Visible','Off');
                set(handles.text30,'Visible','Off');
            elseif handles.popindex==6
                set(handles.pushbutton1,'Visible','On');
                set(handles.pushbutton3,'Visible','On');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','On');
                set(handles.text4,'Visible','On');
                set(handles.text5,'Visible','On');
                set(handles.edit1,'Visible','On');
                set(handles.edit2,'Visible','On');
                set(handles.edit3,'Visible','On');
                set(handles.edit5,'Visible','Off');
                set(handles.text7,'Visible','On');
                set(handles.text8,'Visible','On');
                set(handles.text11,'Visible','Off');
                set(handles.text12,'Visible','Off');
                set(handles.text13,'Visible','Off');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','Off');
                set(handles.edit7,'Visible','Off');
                set(handles.edit8,'Visible','Off');
                set(handles.edit9,'Visible','Off');
                set(handles.text19,'Visible','On');
                set(handles.text20,'Visible','On');
                set(handles.text23,'Visible','Off');
                set(handles.text24,'Visible','On');
                set(handles.text25,'Visible','On');
                set(handles.text26,'Visible','On');
                set(handles.popupmenu4,'Visible','On');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','Off');
                set(handles.edit14,'Visible','Off');
                set(handles.edit15,'Visible','Off');
                set(handles.text28,'Visible','Off');
                set(handles.text29,'Visible','Off');
                set(handles.text30,'Visible','Off');
            elseif handles.popindex==8
                set(handles.pushbutton1,'Visible','On');
                set(handles.pushbutton3,'Visible','On');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','On');
                set(handles.text4,'Visible','On');
                set(handles.text5,'Visible','On');
                set(handles.text7,'Visible','On');
                set(handles.text8,'Visible','On');
                set(handles.edit1,'Visible','On');
                set(handles.edit2,'Visible','On');
                set(handles.edit3,'Visible','On');
                set(handles.edit5,'Visible','On');
                set(handles.text11,'Visible','On');
                set(handles.text12,'Visible','On');
                set(handles.text13,'Visible','On');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','On');
                set(handles.edit7,'Visible','Off');
                set(handles.edit8,'Visible','Off');
                set(handles.edit9,'Visible','Off');
                set(handles.text19,'Visible','On');
                set(handles.text20,'Visible','On');
                set(handles.text23,'Visible','On');
                set(handles.text24,'Visible','On');
                set(handles.text25,'Visible','On');
                set(handles.text26,'Visible','On');
                set(handles.popupmenu4,'Visible','On');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','Off');
                set(handles.edit14,'Visible','Off');
                set(handles.edit15,'Visible','Off');
                set(handles.text28,'Visible','Off');
                set(handles.text29,'Visible','Off');
                set(handles.text30,'Visible','Off');
                set(app.WeightDropDown,'Visible','On');
                set(app.WeightDropDownLabel,'Visible','On')
                set(app.WeightSettingEditField,'Visible','On')
                set(app.WeightSettingEditFieldLabel,'Visible','On')
            elseif handles.popindex==9
                set(handles.pushbutton1,'Visible','On');
                set(handles.pushbutton3,'Visible','On');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','On');
                set(handles.text4,'Visible','On');
                set(handles.text5,'Visible','On');
                set(handles.edit1,'Visible','On');
                set(handles.edit2,'Visible','On');
                set(handles.edit3,'Visible','On');
                set(handles.edit5,'Visible','Off');
                set(handles.text7,'Visible','On');
                set(handles.text8,'Visible','On');
                set(handles.text11,'Visible','On');
                set(handles.text12,'Visible','On');
                set(handles.text13,'Visible','Off');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','On');
                set(handles.edit7,'Visible','Off');
                set(handles.edit8,'Visible','Off');
                set(handles.edit9,'Visible','Off');
                set(handles.text19,'Visible','On');
                set(handles.text20,'Visible','On');
                set(handles.text23,'Visible','On');
                set(handles.text24,'Visible','On');
                set(handles.text25,'Visible','On');
                set(handles.text26,'Visible','On');
                set(handles.popupmenu4,'Visible','On');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','Off');
                set(handles.edit14,'Visible','Off');
                set(handles.edit15,'Visible','Off');
                set(handles.text28,'Visible','Off');
                set(handles.text29,'Visible','Off');
                set(handles.text30,'Visible','Off');
                set(app.WeightDropDown,'Visible','On');
                set(app.WeightDropDownLabel,'Visible','On')
                set(app.WeightSettingEditField,'Visible','On')
                set(app.WeightSettingEditFieldLabel,'Visible','On')
            elseif handles.popindex==7
                set(handles.pushbutton1,'Visible','On');
                set(handles.pushbutton3,'Visible','On');
                set(handles.pushbutton4,'Visible','Off');
                set(handles.text3,'Visible','On');
                set(handles.text4,'Visible','Off');
                set(handles.text5,'Visible','On');
                set(handles.edit1,'Visible','On');
                set(handles.edit2,'Visible','Off');
                set(handles.edit3,'Visible','On');
                set(handles.edit5,'Visible','Off');
                set(handles.text7,'Visible','On');
                set(handles.text8,'Visible','On');
                set(handles.text11,'Visible','On');
                set(handles.text12,'Visible','On');
                set(handles.text13,'Visible','Off');
                set(handles.text15,'Visible','Off');
                set(handles.edit6,'Visible','Off');
                set(handles.text16,'Visible','Off');
                set(handles.edit7,'Visible','Off');
                set(handles.edit8,'Visible','Off');
                set(handles.edit9,'Visible','Off');
                set(handles.text19,'Visible','On');
                set(handles.text20,'Visible','On');
                set(handles.text23,'Visible','Off');
                set(handles.text24,'Visible','On');
                set(handles.text25,'Visible','On');
                set(handles.text26,'Visible','On');
                set(handles.popupmenu4,'Visible','On');
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
                set(handles.pushbutton6,'Visible','On');
                set(handles.edit13,'Visible','On');
                set(handles.edit14,'Visible','On');
                set(handles.edit15,'Visible','On');
                set(handles.text28,'Visible','On');
                set(handles.text29,'Visible','On');
                set(handles.text30,'Visible','On');
            end
            guidata(hObject, handles);
        end

        % Value changed function: popupmenu4
        function popupmenu4_Callback(app, event)
            % --- Executes on selection change in popupmenu4.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to popupmenu4 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles.intindex=get(handles.popupmenu4,'Value');
            p=get(handles.popupmenu4,'Value');
            if p==2
                set(handles.popupmenu5,'Visible','On');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','On');
            elseif p==3
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','On');
                set(handles.text27,'Visible','On');
            elseif p==1
                set(handles.popupmenu5,'Visible','Off');
                set(handles.popupmenu6,'Visible','Off');
                set(handles.text27,'Visible','Off');
            end
            guidata(hObject, handles);
        end

        % Value changed function: popupmenu5
        function popupmenu5_Callback(app, event)
            % --- Executes on selection change in popupmenu5.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to popupmenu5 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            p=get(handles.popupmenu5,'Value');
            if p==2
                handles.accindex=1e-8;
            elseif p==3
                handles.accindex=1e-10;
            elseif p==4
                handles.accindex=1e-12;
            elseif p==5
                handles.accindex=1e-14;
            end
            guidata(hObject, handles);
        end

        % Value changed function: popupmenu6
        function popupmenu6_Callback(app, event)
            % --- Executes on selection change in popupmenu6.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to popupmenu6 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            p=get(handles.popupmenu6,'Value');
            if p==2
                handles.accindex=5e6;
            elseif p==3
                handles.accindex=1e7;
            elseif p==4
                handles.accindex=5e7;
            elseif p==5
                handles.accindex=1e8;
            end
            guidata(hObject, handles);
        end

        % Value changed function: popupmenu7
        function popupmenu7_Callback(app, event)
            % --- Executes on selection change in popupmenu7.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to popupmenu7 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            index2=get(handles.popupmenu7,'Value');
            handles.s=handles.records{1,index2};
            set(handles.edit4,'String',handles.s);
            guidata(hObject, handles);
        end

        % Button pushed function: pushbutton1
        function pushbutton1_Callback(app, event)
            % --- Executes on button press in pushbutton1.
            %SOG---MAIN
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton1 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles.popindex=get(handles.popupmenu1,'Value');
            if handles.popindex==2 || handles.popindex==6
                acc=handles.accindex;
                x=str2double(get(handles.edit3,'string'));handles.xmax=x;
                n=str2double(get(handles.edit1,'string'));handles.n=n;
                T=str2double(get(handles.edit16,'string'));handles.T=T;
                if handles.popindex==2
                    v=eval(get(handles.edit8,'string'));handles.v=v;
                    stepsize=eval(get(handles.edit9,'string'));handles.stepsize=stepsize;
                else
                    v=0;stepsize=0.002;
                    handles.xmin=v;handles.stepsize=stepsize;
                end
                scale=str2double(get(handles.edit2,'string'));handles.scale=scale;
                s=get(handles.edit4,'string');handles.s=s;handles.records=[handles.records,{s}];
                set(handles.popupmenu7,'String',handles.records);
                set(handles.text8,'string','Now Loading');
                if handles.radioindex==0
                    handles.err=exptexting2(app, x,n,scale,s,v,stepsize,acc);
                else
                    s1=s;s2=get(handles.edit10,'string');
                    sp=eval(get(handles.edit11,'string'));
                    handles.err=exptexting5(app, x,n,scale,s1,s2,sp,v,stepsize,acc);
                end
                digitts=max(handles.err(3,:));
                set(handles.text8,'string',num2str(digitts));
                axes(handles.axes1);
                plot(handles.err(1,:),handles.err(3,:));
                xlabel('x');ylabel('Error');
                guidata(hObject, handles);
            elseif handles.popindex==4 || handles.popindex==5 || handles.popindex==8 || handles.popindex==9
                u=eval(get(handles.edit3,'string'));handles.xmax=u;
                n=eval(get(handles.edit1,'string'));handles.n=n;
                T=str2double(get(handles.edit16,'string'));handles.T=T;
                zz=eval(get(handles.edit2,'string'));handles.scale=zz;
                s=get(handles.edit4,'string');handles.s=[s,'--SOG'];handles.records=[handles.records,{s}];
                set(handles.popupmenu7,'String',handles.records);
                g=@(x) eval(s);
                if handles.popindex==4 || handles.popindex==5
                    digitts=eval(get(handles.edit7,'string'));
                    handles.digitss=digitts;
                else
                    digitts=300;
                end
                mp.Digits(digitts);
                scale=n*zz;
                set(handles.text12,'string','Computing Integrals');drawnow;
                if handles.radioindex==0
                    A=zeros(1,2*n);
                    if handles.intindex==2
                        for i=0:2*n-1
                            t=i;
                            set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                            q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
                            A(i+1)=quadgk(q,0,pi,'AbsTol',handles.accindex,'MaxIntervalCount',1e10,'RelTol',0);
                        end
                    elseif handles.intindex==3
                        w=handles.accindex;
                        ww=ceil(sqrt(w));
                        d=2;
                        S=2*ww/d;
                        N=100;
                        [b1,b2]=grule(app, N);
                        for i=0:2*n-1
                            t=i;
                            set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                            q=@(x) (2/pi)*cos(t*x).*g(sqrt(-scale*log((1+cos(x))/2)));
                            for tt=1:S
                                A(i+1)=A(i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
                            end
                            A(i+1)=A(i+1)/ww;
                        end
                    end
                else
                    s1=s;s2=get(handles.edit10,'string');
                    sp=eval(get(handles.edit11,'string'));
                    g1=@(x) eval(s1);g2=@(x) eval(s2);sps=acos(2*exp(-sp*sp/scale)-1);
                    A=zeros(1,2*n);
                    if handles.intindex==2
                        for i=0:2*n-1
                            t=i;
                            set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                            q1=@(x) (2/pi)*cos(t*x).*g1(sqrt(-scale*log((1+cos(x))/2)));
                            q2=@(x) (2/pi)*cos(t*x).*g2(sqrt(-scale*log((1+cos(x))/2)));
                            A(i+1)=quadgk(q1,0,sps,'AbsTol',handles.accindex,'MaxIntervalCount',1e10,'RelTol',0)+quadgk(q2,sps,pi,'AbsTol',1.0e-12,'MaxIntervalCount',1e10,'RelTol',0);
                        end
                    elseif handles.intindex==3
                        w=handles.accindex;
                        ww=ceil(sqrt(w));
                        d=2;
                        S=2*ww/d;
                        N=100;
                        [b1,b2]=grule(app, N);
                        for i=0:2*n-1
                            t=i;
                            set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                            q1=@(x) (2/pi)*cos(t*x).*g1(sqrt(-scale*log((1+cos(x))/2)));
                            q2=@(x) (2/pi)*cos(t*x).*g2(sqrt(-scale*log((1+cos(x))/2)));
                            for tt=1:S
                                A(i+1)=A(i+1)+q1(((2*tt-1+b1-ww)/ww)*sps/2+sps/2)*b2'*sps/2;
                                A(i+1)=A(i+1)+q2(((2*tt-1+b1-ww)/ww)*(pi-sps)/2+(sps+pi)/2)*b2'*(pi-sps)/2;
                            end
                            A(i+1)=A(i+1)/ww;
                        end
                    end
                end
                A=mp(A,digitts);
                maxC=0.00;
                set(handles.text12,'string','Computing Linear Coefficients(VP)');drawnow;
                for r=0:2*n-1
                    set(handles.text11,'string',[num2str(r/2/n*100),'%(2/12)']);drawnow;
                    if(r==0)
                        sum1=mp('0.00');
                        for k=1:n
                            sum1=mp(sum1+mp(-1)^mp(k)*mp(A(k+1)),digitts);
                        end
                        sum2=mp('0.00');
                        for k=1:n-1
                            sum2=mp(sum2+mp(-1)^mp(n+k)*(1-mp(k)/n)*mp(A(n+k+1)),digitts);
                        end
                        C(r+1)=mp(A(r+1)/2+sum1+sum2,digitts);
                    end
                    if(1<=r&&r<=n)
                        sum1=mp('0.00');
                        sum2=mp('0.00');
                        for k=r:n
                            sum1=mp(sum1+mp((-1)^(k-r))*mp(2*mp(k)/mp(k+r))* mp(factorial(mp(k+r))) / mp(factorial(mp(2*r))*factorial(mp(k-r)))*mp((mp(2)^mp(2*r-1)))*mp(A(k+1)),digitts );
                        end
                        for k=1:n-1
                            sum2=mp(sum2+mp((-1)^(n+k-r))*mp(1-mp(k)/n)*2*mp((n+k)*(mp(1)/(n+k+r)))*(mp(factorial(mp(n+k+r))))/((mp(factorial(mp(2*r))))*(mp(factorial(mp(n+k-r)))))*mp(2^mp(2*r-1))*mp(A(n+k+1)),digitts);
                        end
                        C(r+1)=mp((sum1)+(sum2),digitts);
                    end
                    if(r>=n+1)
                        sum1=mp('0.00');
                        for k=(r-n:n-1)
                            sum1=mp(sum1+mp(mp(-1)^mp(n+k-r))*mp(1-mp(k)/n)*2*mp(mp(n+k)/(mp(n+k+r)))*(mp(factorial(mp(n+k+r))))/(mp(factorial(mp(2*r)))*mp(factorial(mp(n+k-r))))*(mp(2^mp(2*r-1)))*mp(A(n+k+1)),digitts);
                        end
                        C(r+1)=mp(sum1,digitts);
                    end
                end
                if handles.popindex==4 || handles.popindex==5
                    v=eval(get(handles.edit8,'string'));handles.xmin=v;
                    stepsize=eval(get(handles.edit9,'string'));handles.stepsize=stepsize;
                    x=(v+stepsize):stepsize:u;
                else
                    stepsize=0.002;v=0;
                    handles.xmin=v;handles.stepsize=stepsize;
                    if u>=10
                        stepsize=0.2;
                    end
                    x=stepsize:stepsize:u;
                end
                if handles.radioindex2==0
                    NNN=size(x,2);
            %         for i=1:NNN
            %             i;
            %             a(i)=mp('0.00');
            %         end
                    set(handles.text12,'string','Testing the VP sum');drawnow;
                    set(handles.text11,'string','(3/12)');drawnow;
            %         for k=0:2*n-1
            %             set(handles.text11,'string',[num2str(k/2/n*100),'%(3/12)']);drawnow;
            %             for i=1:NNN
            %                 y(i)=mp(exp(mp(-mp(k)*(mp(x(i))*mp(x(i)))/mp(scale))));
            %                 a(i)=mp(a(i)+mp(C(k+1))*mp(y(i)),digitts);
            %             end
            %         end
                    CC=C.';xx=mp(x).^2;KK=mp(0:(2*n-1))/mp(scale);
                    a=sum(CC.*exp(-KK'.*xx),1);
                    eta=double(max(mp(abs(mp(a-(g(x)))))));
                    if handles.radioindex==1
                        eta1=[ones(1,round((sp-v)/stepsize)) zeros(1,round(NNN-(sp-v)/stepsize))];
                        eta2=[zeros(1,round((sp-v)/stepsize)) ones(1,round(NNN-(sp-v)/stepsize))];
                        gg=g1(x).*eta1+g2(x).*eta2;
                        eta=double(max(mp(abs(mp(a-gg)))));
                    end
                    set(handles.text8,'string',[num2str(eta) '(VP)']);drawnow;
                end
                X=C;handles.XX=double(X(1));
                %%Model Reduction
                set(handles.text12,'string','Model Reduction');drawnow;
                set(handles.text11,'string','(4/12)');drawnow;
                n=2*n;mp.Digits(digitts);
                if handles.popindex==4 || handles.popindex==5
                    v=eval(get(handles.edit8,'string'));
                    stepsize=eval(get(handles.edit9,'string'));
                    x=(v+stepsize):stepsize:u;
                else
                    stepsize=0.002;
                    if u>=10
                        stepsize=0.2;
                    end
                    x=stepsize:stepsize:u;
                end
                N=size(x,2);
                y=g(x);
                if handles.radioindex==1
                    NNN=size(x,2);
                    eta1=[ones(1,round((sp-v)/stepsize)) zeros(1,round(NNN-(sp-v)/stepsize))];
                    eta2=[zeros(1,round((sp-v)/stepsize)) ones(1,round(NNN-(sp-v)/stepsize))];
                    y=g1(x).*eta1+g2(x).*eta2;
                end
                y=y.';
                A=mp(-mp(diag(1:n-1))/mp(n*zz/2),digitts);
                A2=mp(-mp((1:n-1))/mp(n*zz/2),digitts);
                EA2 = mp(exp(A2*T),digitts);
                B=mp(zeros(n-1,1),digitts);
                C=mp(zeros(1,n-1),digitts);
                for i=1:n-1
                    if(X(i+1)>0)
                        C(i)=mp(sqrt(mp(X(i+1),digitts)),digitts);
                        B(i)=mp(sqrt(mp(X(i+1),digitts)),digitts);
                    else
                        C(i)=mp(sqrt(mp(-X(i+1),digitts)),digitts);
                        B(i)=mp(-sqrt(mp(-X(i+1),digitts)),digitts);
                    end
                end
                BBB=B*B.';
                CCC=C.'*C;
                AA=A2+A2.';
                set(handles.text12,'string','Preparing P and Q');drawnow;
                set(handles.text11,'string',' ');drawnow;
                p=app.WeightDropDown.Value;
                handles.weight ='Empty';
                handles.weightset = [];
                if strcmp(p, '1')
                    handles.weight ='1';
                    handles.weightset = [];
                    EAA = mp(mp(ones(n-1,n-1),digitts) - mp(EA2.'*EA2,digitts),digitts);
                    P=-BBB.*(EAA./AA);Q=-CCC.*(EAA./AA);
                elseif strcmp(p, '(x+a).^(-1/2)')
                    a=eval(get(app.WeightSettingEditField,'Value'));
                    handles.weight ='(x+a).^(-1/2)';
                    handles.weightset = string(a);
                    a = mp(a,digitts);
                    T = mp(T,digitts);
                    EAA = mp(expint(mp(-a.*AA,digitts)),digitts) - mp(expint(mp(-(T+a).*AA,digitts)),digitts);
                    EAA = mp(mp(exp(mp(a.*AA,digitts)),digitts).*EAA,digitts);
                    P=-BBB.*EAA;Q=-CCC.*EAA;
                elseif strcmp(p, 'exp(-s.*x)')
                    a=eval(get(app.WeightSettingEditField,'Value'));
                    handles.weight ='exp(-s.*x)';
                    handles.weightset = string(a);
                    a = mp(a,digitts);
                    T = mp(T,digitts);
                    AA = AA - 2*a;
                    EA2 = mp(exp((A2-a)*T),digitts);
                    EAA = mp(ones(n-1,n-1),digitts) - EA2.'*EA2;
                    P=-BBB.*(EAA./AA);Q=-CCC.*(EAA./AA);
                elseif strcmp(p, 'Custom')
                    wf = get(app.WeightSettingEditField,'Value');
                    handles.weight ='Custom';
                    handles.weightset = wf;
                    weight_func = @(x)eval(wf);
                    Int_func = @(x)weight_func(x).^2.*exp(AA.*x);
                    EAA = integral(Int_func,0,T,"ArrayValued",true,"AbsTol",1e-16,"RelTol",0);
                    P=-BBB.*EAA;Q=-CCC.*EAA;
                end
                %P(1,1)=-BBB(1,1)*T;Q(1,1)=-CCC(1,1)*T;
                set(handles.text12,'string','Cholesky');drawnow;
                mp.Digits(digitts);
                for i=1:n-1
                    set(handles.text11,'string',[num2str(i/n*100),'%(6/12)']);drawnow;
                    P(i,i)=sqrt(P(i,i)-sum(P(i,1:i-1).^2));
                    Q(i,i)=sqrt(Q(i,i)-sum(Q(i,1:i-1).^2));
                    PP=sum(P(i+1:n-1,1:i-1).*P(i,1:i-1),2);
                    QQ=sum(Q(i+1:n-1,1:i-1).*Q(i,1:i-1),2);
                    P(i+1:n-1,i)=(P(i+1:n-1,i)-PP(1:n-1-i))/P(i,i);
                    Q(i+1:n-1,i)=(Q(i+1:n-1,i)-QQ(1:n-1-i))/Q(i,i);
                end
                P=tril(ones(n-1)).*P;Q=tril(ones(n-1)).*Q;
                Lc=P;
                LL=P.'*Q;
                set(handles.text12,'string','SVD.It takes some time.');drawnow;
                set(handles.text11,'string','(7/12)');drawnow;
                [U,sigma,V]=svd(LL);
                set(handles.text12,'string','Preparing LLL');drawnow;
                set(handles.text11,'string','(8/12)');drawnow;
                LLL=mp(diag(mp(sigma,digitts)),digitts);
                LLL=mp(LLL.^(mp(-1/2,digitts)),digitts);
                set(handles.text12,'string','Preparing T');drawnow;
                set(handles.text11,'string','(9/12)');drawnow;
                T=mp(Lc,digitts)*mp(U,digitts)*mp(diag(LLL),digitts);
                set(handles.text12,'string','Preparing Anew');drawnow;
                set(handles.text11,'string','(10/12)');drawnow;
                Anew=mp(inv(mp(T,digitts)),digitts)*mp(A,digitts)*mp(T,digitts);
                set(handles.text12,'string','Preparing Bnew');drawnow;
                set(handles.text11,'string','(11/12)');drawnow;
                Bnew=mp(inv(mp(T,digitts)),digitts)*mp(B,digitts);
                Cnew=mp(C,digitts)*mp(T,digitts);
                handles.Anew=Anew;
                handles.Bnew=Bnew;
                handles.Cnew=Cnew;
                if handles.popindex==4 || handles.popindex==8
                    choos=str2double(get(handles.edit5,'string'));
                    handles.mrterms=choos;
                    Anew=Anew(1:choos,1:choos);
                    Bnew=Bnew(1:choos);
                    Cnew=Cnew(1:choos);
                    Anew_double=double(Anew);
                    [Avec Aeig]=eig(Anew_double);
                    Bnew_double=pinv(Avec,1e-14)*Bnew;
                    Cnew_double=Cnew*Avec;
                    for i=1:choos
                        w(i)=Bnew_double(i)*Cnew_double(i);
                    end
                    p=diag(Aeig);
                    approx=zeros(N,1);
                    approx=approx + double(X(1));
                    if handles.popindex==4 || handles.popindex==5
                        v=eval(get(handles.edit8,'string'));
                        stepsize=eval(get(handles.edit9,'string'));
                        x=(v+stepsize):stepsize:u;
                    else
                        stepsize=0.002;
                        if u>=10
                            stepsize=0.2;
                        end
                        x=stepsize:stepsize:u;
                    end
                    set(handles.text12,'string','Testing MR');drawnow;
                    for i=1:N
                        set(handles.text11,'string',[num2str(i/N*100),'%(12/12)']);drawnow;
                        for j=1:choos
                            approx(i)=approx(i)+w(j)*exp((p(j))*x(i)^2);
                        end
                    end
                    approx=double(approx);
                    errors=abs(approx-y);
                    handles.errors=max(abs(errors));
                    axes(handles.axes1);
                    plot(x,errors);
                    xlabel('x');ylabel('Error');
                    set(handles.text8,'string',num2str(max(abs(approx-y))));drawnow;
                    handles.pp=double(p);
                    handles.ww=double(w);
                    set(handles.text12,'string','Clear!Please save.');drawnow;
                    set(handles.text11,'string','100%');drawnow;
                    set(handles.pushbutton4,'Visible','On');
                    set(handles.text15,'Visible','On');
                    set(handles.edit6,'Visible','On');
                    guidata(hObject, handles);
                else
                    if handles.popindex==5
                        v=eval(get(handles.edit8,'string'));
                        stepsize=eval(get(handles.edit9,'string'));
                        x=(v+stepsize):stepsize:u;
                    else
                        stepsize=0.002;
                        if u>=10
                            stepsize=0.2;
                        end
                        x=stepsize:stepsize:u;
                    end
                    errs=zeros(1,n);
                    xx=x.^2;
                    for choos=1:n-1
                        set(handles.text12,'string','Testing MR terms');drawnow;
                        set(handles.text11,'string',[num2str(choos/n*100) '%(11/11)']);drawnow;
                        Anews=Anew(1:choos,1:choos);
                        Bnews=Bnew(1:choos);
                        Cnews=Cnew(1:choos);
                        Anew_double=double(Anews);
                        [Avec Aeig]=eig(Anew_double);
                        Bnew_double=pinv(Avec,1e-14)*Bnews;
                        Cnew_double=Cnews*Avec;
                        w=zeros(choos,1);
                        for i=1:choos
                            w(i)=Bnew_double(i)*Cnew_double(i);
                        end
                        p=diag(Aeig);
            %             approx=zeros(N,1);
            %             approx=approx+double(X(1));
            %             for i=1:N
            %                 for j=1:choos
            %                     approx(i)=approx(i)+w(j)*exp((p(j))*x(i)^2);
            %                 end
            %             end
                        approx=sum(w.*exp(p.*xx),1)+double(X(1));
            %             errs(choos)=max(abs(approx-y));
                        errs(choos)=max(abs(approx-y.'));
                    end
                    handles.errs=errs;
                    axes(handles.axes1);
                    semilogy(1:n,errs);
                    xlabel('The number of terms');ylabel('Maximum error');
                    set(handles.text12,'string','Clear!');drawnow;
                    set(handles.text11,'string','100%');drawnow;
                    set(handles.text13,'Visible','On');
                    set(handles.edit5,'Visible','On');
                    set(handles.pushbutton4,'Visible','On');
                    set(handles.text15,'Visible','On');
                    set(handles.edit6,'Visible','On');
                    guidata(hObject, handles);
                end
            elseif handles.popindex==3 || handles.popindex==7
                acc=handles.accindex;
                x=str2double(get(handles.edit3,'string'));handles.xmax=x;
                n=str2double(get(handles.edit1,'string'));handles.n=n;
                T=eval(get(handles.edit16,'string'));handles.T=T;
                if handles.popindex==3
                    v=eval(get(handles.edit8,'string'));handles.v=v;
                    stepsize=eval(get(handles.edit9,'string'));handles.stepsize=stepsize;
                else
                    v=0;stepsize=0.002;
                    handles.xmin=v;handles.stepsize=stepsize;
                end
                s=get(handles.edit4,'string');handles.s=s;handles.records=[handles.records,{s}];
                set(handles.popupmenu7,'String',handles.records);
                scstart=eval(get(handles.edit13,'string'));
                scend=eval(get(handles.edit14,'string'));
                scstep=eval(get(handles.edit15,'string'));
                scale=scstart:scstep:scend;
                handles.err=mexptexting2(app, x,n,scale,s,v,stepsize,acc);
                axes(handles.axes1);
                semilogy(scale,handles.err);
                xlabel('scale/n');ylabel('Error');
                err=min(handles.err);
                whereerr=find(handles.err==err);
                set(handles.text8,'string',num2str(err));drawnow;
                set(handles.text11,'string',['scale/n=' num2str(scale(whereerr(1)))]);drawnow;
                guidata(hObject, handles);
            end
        end

        % Button pushed function: pushbutton3
        function pushbutton3_Callback(app, event)
            % --- Executes on button press in pushbutton3.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton3 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles.popindex=get(handles.popupmenu1,'Value');
            if handles.popindex==2 || handles.popindex==6
                acc=handles.accindex;
                x=str2double(get(handles.edit3,'string'));
                n=str2double(get(handles.edit1,'string'));
                T = str2double(app.edit16.Value);handles.T=T;
                if handles.popindex==2
                    v=eval(get(handles.edit8,'string'));
                    stepsize=eval(get(handles.edit9,'string'));
                else
                    v=0;stepsize=0.002;
                end
                scale=str2double(get(handles.edit2,'string'));
                s=get(handles.edit4,'string');handles.records=[handles.records,{s}];handles.s=s;
                set(handles.popupmenu7,'String',handles.records);
                set(handles.text8,'string','Now Loading');
                if handles.radioindex==0
                    handles.err=exptexting3(app, x,n,scale,s,v,stepsize,acc);
                else
                    s1=s;s2=get(handles.edit10,'string');
                    sp=eval(get(handles.edit11,'string'));
                    handles.err=exptexting6(app, x,n,scale,s1,s2,sp,v,stepsize,acc);
                end
                digitts=max(handles.err(3,:));
                set(handles.text8,'string',num2str(digitts));
                axes(handles.axes1);
                plot(handles.err(1,:),handles.err(3,:));
                xlabel('x');ylabel('Error');
                guidata(hObject, handles);
            elseif handles.popindex==4 || handles.popindex==5 || handles.popindex==8 || handles.popindex==9
                u=eval(get(handles.edit3,'string'));handles.xmax=u;
                n=eval(get(handles.edit1,'string'));handles.n=n;
                T=str2double(app.edit16.Value);handles.T=T;
                zz=eval(get(handles.edit2,'string'));handles.scale=zz;
                s=get(handles.edit4,'string');handles.s=[s,'--SOE'];handles.records=[handles.records,{s}];
                set(handles.popupmenu7,'String',handles.records);
                gg=@(x) eval(s);
                g=@(x)truncation(app,x,gg);
                if handles.popindex==4 || handles.popindex==5
                    digitts=eval(get(handles.edit7,'string'));
                    handles.digitss=digitts;
                else
                    digitts=300;
                end
                mp.Digits(digitts);
                scale=n*zz;
                set(handles.text12,'string','Computing Integrals');drawnow;
                if handles.radioindex==0
                    A=zeros(1,2*n);
                    if handles.intindex==2
                        for i=0:2*n-1
                            t=i;
                            set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                            q=@(x) (2/pi)*cos(t*x).*g(-scale*log((1+cos(x))/2));
                            A(i+1)=quadgk(q,0,pi,'AbsTol',handles.accindex,'MaxIntervalCount',1e10,'RelTol',0);
                        end
                    elseif handles.intindex==3
                        w=handles.accindex;
                        ww=ceil(sqrt(w));
                        d=2;
                        S=2*ww/d;
                        N=100;
                        [b1,b2]=grule(app, N);
                        for i=0:2*n-1
                            t=i;
                            set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                            q=@(x) (2/pi)*cos(t*x).*g(-scale*log((1+cos(x))/2));
                            for tt=1:S
                                A(i+1)=A(i+1)+q(((2*tt-1+b1-ww)/ww)*pi/2+pi/2)*b2'*pi/2;
                            end
                            A(i+1)=A(i+1)/ww;
                        end
                    end
                else
                    s1=s;s2=get(handles.edit10,'string');
                    sp=eval(get(handles.edit11,'string'));
                    g1=@(x) eval(s1);g2=@(x) eval(s2);sps=acos(2*exp(-sp/scale)-1);
                    A=zeros(1,2*n);
                    if handles.intindex==2
                        for i=0:2*n-1
                            t=i;
                            set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                            q1=@(x) (2/pi)*cos(t*x).*g1((-scale*log((1+cos(x))/2)));
                            q2=@(x) (2/pi)*cos(t*x).*g2((-scale*log((1+cos(x))/2)));
                            A(i+1)=quadgk(q1,0,sps,'AbsTol',handles.accindex,'MaxIntervalCount',1e10,'RelTol',1e-10)+quadgk(q2,sps,pi,'AbsTol',1.0e-12,'MaxIntervalCount',1e10,'RelTol',1e-10);
                        end
                    elseif handles.intindex==3
                        w=handles.accindex;
                        ww=ceil(sqrt(w));
                        d=2;
                        S=2*ww/d;
                        N=100;
                        [b1,b2]=grule(app, N);
                        for i=0:2*n-1
                            t=i;
                            set(handles.text11,'string',[num2str(t/2/n*100),'%(1/12)']);drawnow;
                            q1=@(x) (2/pi)*cos(t*x).*g1((-scale*log((1+cos(x))/2)));
                            q2=@(x) (2/pi)*cos(t*x).*g2((-scale*log((1+cos(x))/2)));
                            for tt=1:S
                                A(i+1)=A(i+1)+q1(((2*tt-1+b1-ww)/ww)*sps/2+sps/2)*b2'*sps/2;
                                A(i+1)=A(i+1)+q2(((2*tt-1+b1-ww)/ww)*(pi-sps)/2+(sps+pi)/2)*b2'*(pi-sps)/2;
                            end
                            A(i+1)=A(i+1)/ww;
                        end
                    end
                end
                A=mp(A,digitts);
                maxC=0.00;
                set(handles.text12,'string','Computing Linear Coefficients(VP)');drawnow;
                for r=0:2*n-1
                    set(handles.text11,'string',[num2str(r/2/n*100),'%(2/12)']);drawnow;
                    if(r==0)
                        sum1=mp('0.00');
                        for k=1:n
                            sum1=mp(sum1+mp(-1)^mp(k)*mp(A(k+1)),digitts);
                        end
                        sum2=mp('0.00');
                        for k=1:n-1
                            sum2=mp(sum2+mp(-1)^mp(n+k)*(1-mp(k)/n)*mp(A(n+k+1)),digitts);
                        end
                        C(r+1)=mp(A(r+1)/2+sum1+sum2,digitts);
                    end
                    if(1<=r&&r<=n)
                        sum1=mp('0.00');
                        sum2=mp('0.00');
                        for k=r:n
                            sum1=mp(sum1+mp((-1)^(k-r))*mp(2*mp(k)/mp(k+r))* mp(factorial(mp(k+r))) / mp(factorial(mp(2*r))*factorial(mp(k-r)))*mp((mp(2)^mp(2*r-1)))*mp(A(k+1)),digitts );
                        end
                        for k=1:n-1
                            sum2=mp(sum2+mp((-1)^(n+k-r))*mp(1-mp(k)/n)*2*mp((n+k)*(mp(1)/(n+k+r)))*(mp(factorial(mp(n+k+r))))/((mp(factorial(mp(2*r))))*(mp(factorial(mp(n+k-r)))))*mp(2^mp(2*r-1))*mp(A(n+k+1)),digitts);
                        end
                        C(r+1)=mp((sum1)+(sum2),digitts);
                    end
                    if(r>=n+1)
                        sum1=mp('0.00');
                        for k=(r-n:n-1)
                            sum1=mp(sum1+mp(mp(-1)^mp(n+k-r))*mp(1-mp(k)/n)*2*mp(mp(n+k)/(mp(n+k+r)))*(mp(factorial(mp(n+k+r))))/(mp(factorial(mp(2*r)))*mp(factorial(mp(n+k-r))))*(mp(2^mp(2*r-1)))*mp(A(n+k+1)),digitts);
                        end
                        C(r+1)=mp(sum1,digitts);
                    end
                end
                if handles.popindex==4 || handles.popindex==5
                    v=eval(get(handles.edit8,'string'));
                    stepsize=eval(get(handles.edit9,'string'));
                    handles.xmin=v;handles.stepsize=stepsize;
                    x=(v+stepsize):stepsize:u;
                else
                    stepsize=0.002;v=0;
                    handles.xmin=v;handles.stepsize=stepsize;
                    if u>=10
                        stepsize=0.2;
                    end
                    x=stepsize:stepsize:u;
                end
                if handles.radioindex2==0
                    NNN=size(x,2);
            %         for i=1:NNN
            %             i;
            %             a(i)=mp('0.00');
            %         end
                    set(handles.text12,'string','Testing the VP sum');drawnow;
            %         for k=0:2*n-1
            %             set(handles.text11,'string',[num2str(k/2/n*100),'%(2/11)']);drawnow;
            %             for i=1:NNN
            %                 y(i)=mp(exp(mp(-mp(k)*(mp(x(i)))/mp(scale))));
            %                 a(i)=mp(a(i)+mp(C(k+1))*mp(y(i)),digitts);
            %             end
            %         end
                    set(handles.text11,'string','(3/12)');drawnow;
                    CC=C.';xx=mp(x);KK=mp(0:(2*n-1))/mp(scale);
                    a=sum(CC.*exp(-KK'.*xx),1);
                    eta=double(max(mp(abs(mp(a-(g(x)))))));
            %         eta=double(max(mp(abs(mp(a-(g(x)))))));
                    if handles.radioindex==1
                        eta1=[ones(1,round((sp-v)/stepsize)) zeros(1,round(NNN-(sp-v)/stepsize))];
                        eta2=[zeros(1,round((sp-v)/stepsize)) ones(1,round(NNN-(sp-v)/stepsize))];
                        gg=g1(x).*eta1+g2(x).*eta2;
                        eta=double(max(mp(abs(mp(a-gg)))));
                    end
                    set(handles.text8,'string',[num2str(eta) '(VP)']);drawnow;
                end
                X=C;handles.XX=double(X(1));
                %%Model Reduction
                set(handles.text12,'string','Model Reduction');drawnow;
                set(handles.text11,'string','(4/12)');drawnow;
                n=2*n;mp.Digits(digitts);
                if handles.popindex==4 || handles.popindex==5
                    v=eval(get(handles.edit8,'string'));
                    stepsize=eval(get(handles.edit9,'string'));
                    x=(v+stepsize):stepsize:u;
                else
                    stepsize=0.002;
                    if u>=10
                        stepsize=0.2;
                    end
                    x=stepsize:stepsize:u;
                end
                N=size(x,2);
                y=g(x);
                if handles.radioindex==1
                    y=g1(x).*eta1+g2(x).*eta2;
                end
                y=y.';
                A=mp(-mp(diag(1:n-1))/mp(n*zz/2),digitts);
                A2=mp(-mp((1:n-1))/mp(n*zz/2),digitts);
                B=mp(zeros(n-1,1),digitts);
                C=mp(zeros(1,n-1),digitts);
                for i=1:n-1
                    if(X(i+1)>0)
                        C(i)=mp(sqrt(mp(X(i+1),digitts)),digitts);
                        B(i)=mp(sqrt(mp(X(i+1),digitts)),digitts);
                    else
                        C(i)=mp(sqrt(mp(-X(i+1),digitts)),digitts);
                        B(i)=mp(-sqrt(mp(-X(i+1),digitts)),digitts);
                    end
                end
                BBB=B*B.';
                CCC=C.'*C;
                EA2 = mp(exp(A2*T),digitts);
                AA=mp(A2.'+A2,digitts);
                set(handles.text12,'string','Preparing P and Q');drawnow;
                set(handles.text11,'string',' ');drawnow;
                p=app.WeightDropDown.Value;
                handles.weight ='Empty';
                handles.weightset = [];
                if strcmp(p, '1')
                    handles.weight ='1';
                    handles.weightset = [];
                    EAA = mp(mp(ones(n-1,n-1),digitts) - mp(EA2.'*EA2,digitts),digitts);
                    P=-BBB.*(EAA./AA);Q=-CCC.*(EAA./AA);
                elseif strcmp(p, '(x+a).^(-1/2)')
                    a=eval(get(app.WeightSettingEditField,'Value'));
                    handles.weight ='(x+a).^(-1/2)';
                    handles.weightset = string(a);
                    a = mp(a,digitts);
                    T = mp(T,digitts);
                    EAA = mp(expint(mp(-a.*AA,digitts)),digitts) - mp(expint(mp(-(T+a).*AA,digitts)),digitts);
                    EAA = mp(mp(exp(mp(a.*AA,digitts)),digitts).*EAA,digitts);
                    P=-BBB.*EAA;Q=-CCC.*EAA;
                elseif strcmp(p, 'exp(-s.*x)')
                    a=eval(get(app.WeightSettingEditField,'Value'));
                    handles.weight ='exp(-s.*x)';
                    handles.weightset = string(a);
                    a = mp(a,digitts);
                    T = mp(T,digitts);
                    AA = AA - 2*a;
                    EA2 = mp(exp((A2-a)*T),digitts);
                    EAA = mp(ones(n-1,n-1),digitts) - EA2.'*EA2;
                    P=-BBB.*(EAA./AA);Q=-CCC.*(EAA./AA);
                elseif strcmp(p, 'Custom')
                    wf = get(app.WeightSettingEditField,'Value');
                    handles.weight ='Custom';
                    handles.weightset = wf;
                    weight_func = @(x)eval(wf);
                    Int_func = @(x)weight_func(x).^2.*exp(AA.*x);
                    EAA = integral(Int_func,0,T,"ArrayValued",true,"AbsTol",1e-16,"RelTol",0);
                    P=-BBB.*EAA;Q=-CCC.*EAA;
                end
                %P(1,1)=-BBB(1,1)*T;Q(1,1)=-CCC(1,1)*T;
                set(handles.text12,'string','Cholesky');drawnow;
                mp.Digits(digitts);
                for i=1:n-1
                    set(handles.text11,'string',[num2str(i/n*100),'%(6/12)']);drawnow;
                    P(i,i)=sqrt(P(i,i)-sum(P(i,1:i-1).^2));
                    Q(i,i)=sqrt(Q(i,i)-sum(Q(i,1:i-1).^2));
                    PP=sum(P(i+1:n-1,1:i-1).*P(i,1:i-1),2);
                    QQ=sum(Q(i+1:n-1,1:i-1).*Q(i,1:i-1),2);
                    P(i+1:n-1,i)=(P(i+1:n-1,i)-PP(1:n-1-i))/P(i,i);
                    Q(i+1:n-1,i)=(Q(i+1:n-1,i)-QQ(1:n-1-i))/Q(i,i);
                end
                P=tril(ones(n-1)).*P;Q=tril(ones(n-1)).*Q;
                Lc=P;
                LL=P.'*Q;
                set(handles.text12,'string','SVD.It takes some time.');drawnow;
                set(handles.text11,'string','(7/12)');drawnow;
                [U,sigma,V]=svd(LL);
                set(handles.text12,'string','Preparing LLL');drawnow;
                set(handles.text11,'string','(8/12)');drawnow;
                LLL=mp(diag(mp(sigma,digitts)),digitts);
                LLL=mp(LLL.^(mp(-1/2,digitts)),digitts);
                set(handles.text12,'string','Preparing T');drawnow;
                set(handles.text11,'string','(9/12)');drawnow;
                T=mp(Lc,digitts)*mp(U,digitts)*mp(diag(LLL),digitts);
                InvT = mp(diag(LLL),digitts)*mp(V',digitts)*mp(Q',digitts);
                set(handles.text12,'string','Preparing Anew');drawnow;
                set(handles.text11,'string','(10/12)');drawnow;
                Anew=mp(InvT,digitts)*mp(A,digitts)*mp(T,digitts);
                set(handles.text12,'string','Preparing Bnew');drawnow;
                set(handles.text11,'string','(11/12)');drawnow;
                Bnew=mp(InvT,digitts)*mp(B,digitts);
                Cnew=mp(C,digitts)*mp(T,digitts);
                handles.Anew=Anew;
                handles.Bnew=Bnew;
                handles.Cnew=Cnew;
                if handles.popindex==4 || handles.popindex==8
                    choos=str2double(get(handles.edit5,'string'));
                    handles.mrterms=choos;
                    Anew=Anew(1:choos,1:choos);
                    Bnew=Bnew(1:choos);
                    Cnew=Cnew(1:choos);
                    Anew_double=double(Anew);
                    [Avec Aeig]=eig(Anew_double);
                    Bnew_double=pinv(Avec,1e-14)*Bnew;
                    Cnew_double=Cnew*Avec;
                    for i=1:choos
                        w(i)=Bnew_double(i)*Cnew_double(i);
                    end
                    p=diag(Aeig);
                    approx=zeros(N,1);
                    approx=approx + double(X(1));
                    if handles.popindex==4 || handles.popindex==5
                        v=eval(get(handles.edit8,'string'));
                        stepsize=eval(get(handles.edit9,'string'));
                        x=(v+stepsize):stepsize:u;
                    else
                        stepsize=0.002;
                        if u>=10
                            stepsize=0.2;
                        end
                        x=stepsize:stepsize:u;
                    end
                    set(handles.text12,'string','Testing MR');drawnow;
                    for i=1:N
                        set(handles.text11,'string',[num2str(i/N*100),'%(12/12)']);drawnow;
                        for j=1:choos
                            approx(i)=approx(i)+w(j)*exp((p(j))*x(i));
                        end
                    end
                    approx=double(approx);
                    errors=abs(approx-y);
                    axes(handles.axes1);
                    plot(x,errors);
                    xlabel('x');ylabel('Error');
                    handles.errors=max(abs(errors));
                    set(handles.text8,'string',num2str(max(abs(approx-y))));drawnow;
                    handles.pp=double(p);
                    handles.ww=double(w);
                    set(handles.text12,'string','Clear!Please save.');drawnow;
                    set(handles.text11,'string','100%');drawnow;
                    set(handles.pushbutton4,'Visible','On');
                    set(handles.text15,'Visible','On');
                    set(handles.edit6,'Visible','On');
                    guidata(hObject, handles);
                else
                    if handles.popindex==5
                        v=eval(get(handles.edit8,'string'));
                        stepsize=eval(get(handles.edit9,'string'));
                        x=(v+stepsize):stepsize:u;
                    else
                        stepsize=0.002;
                        if u>=10
                            stepsize=0.2;
                        end
                        x=stepsize:stepsize:u;
                    end
                    errs=zeros(1,n);
                    for choos=1:n-1
                        set(handles.text12,'string','Testing MR terms');drawnow;
                        set(handles.text11,'string',[num2str(choos/n*100) '%(11/11)']);drawnow;
                        Anews=Anew(1:choos,1:choos);
                        Bnews=Bnew(1:choos);
                        Cnews=Cnew(1:choos);
                        Anew_double=double(Anews);
                        [Avec Aeig]=eig(Anew_double);
                        Bnew_double=pinv(Avec,1e-14)*Bnews;
                        Cnew_double=Cnews*Avec;
                        w=zeros(choos,1);
                        for i=1:choos
                            w(i)=Bnew_double(i)*Cnew_double(i);
                        end
                        p=diag(Aeig);
            %             approx=zeros(N,1);
            %             approx=approx+double(X(1));
            %             for i=1:N
            %                 for j=1:choos
            %                     approx(i)=approx(i)+w(j)*exp((p(j))*x(i));
            %                 end
            %             end
                        approx=sum(w.*exp(p.*x),1) + double(X(1));
                        errs(choos)=max(abs(approx-y.'));
                    end
                    handles.errs=errs;
                    axes(handles.axes1);
                    semilogy(1:n,errs);
                    xlabel('The number of terms');ylabel('Maximum Error');
                    set(handles.text12,'string','Clear!');drawnow;
                    set(handles.text11,'string','100%');drawnow;
                    set(handles.text13,'Visible','On');
                    set(handles.edit5,'Visible','On');
                    set(handles.pushbutton4,'Visible','On');
                    set(handles.text15,'Visible','On');
                    set(handles.edit6,'Visible','On');
                    guidata(hObject, handles);
                end
            elseif handles.popindex==3 || handles.popindex==7
                acc=handles.accindex;
                x=str2double(get(handles.edit3,'string'));handles.xmax=x;
                n=str2double(get(handles.edit1,'string'));handles.n=n;
                T = str2double(app.edit16.Value);handles.T=T;
                if handles.popindex==7
                    v=eval(get(handles.edit8,'string'));handles.v=v;
                    stepsize=eval(get(handles.edit9,'string'));handles.stepsize=stepsize;
                else
                    v=0;stepsize=0.002;
                    handles.xmin=v;handles.stepsize=stepsize;
                end
                guidata(hObject, handles);
                s=get(handles.edit4,'string');handles.s=s;handles.records=[handles.records,{s}];
                set(handles.popupmenu7,'String',handles.records);
                scstart=eval(get(handles.edit13,'string'));
                scend=eval(get(handles.edit14,'string'));
                scstep=eval(get(handles.edit15,'string'));
                scale=scstart:scstep:scend;
                handles.err=mexptexting3(app, x,n,scale,s,v,stepsize,acc);
                axes(handles.axes1);
                semilogy(scale,handles.err);
                xlabel('scale/n');ylabel('Error');
                err=min(handles.err);
                whereerr=find(handles.err==err);
                set(handles.text8,'string',num2str(err));drawnow;
                set(handles.text11,'string',['scale/n=' num2str(scale(whereerr(1)))]);drawnow;
            end
        end

        % Button pushed function: pushbutton4
        function pushbutton4_Callback(app, event)
            % --- Executes on button press in pushbutton4.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton4 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            Xs=handles.XX;
            ps=handles.pp;
            ws=handles.ww;
            if handles.popindex==5 || handles.popindex==9
                choos=str2num(get(handles.edit5,'string'));
                handles.mrterms=choos;
                if isscalar(choos)
                    Anew=handles.Anew;
                    Bnew=handles.Bnew;
                    Cnew=handles.Cnew;
                    Anews=Anew(1:choos,1:choos);
                    Bnews=Bnew(1:choos);
                    Cnews=Cnew(1:choos);
                    Anew_double=double(Anews);
                    [Avec Aeig]=eig(Anew_double);
                    Bnew_double=pinv(Avec,1e-14)*Bnews;
                    Cnew_double=Cnews*Avec;
                    for i=1:choos
                        ws(i)=Bnew_double(i)*Cnew_double(i);
                    end
                    ps=diag(Aeig);
                    errs=handles.errs;
                    handles.errors=errs(choos);
                else
                    errs=handles.errs;
                    handles.errors = errs(choos);
                end
            end
            pathname=[get(handles.edit6,'string'),'/parameter_',datestr(now,30)];
            mkdir(pathname);
            if isscalar(choos)
                save([pathname,'/X'],'Xs');
                save([pathname,'/s'],'ps');
                save([pathname,'/w'],'ws');
            else
                save([pathname,'/X'],'Xs');
                Anew=handles.Anew;
                Bnew=handles.Bnew;
                Cnew=handles.Cnew;
                for i=choos
                    Anews=Anew(1:i,1:i);
                    Bnews=Bnew(1:i);
                    Cnews=Cnew(1:i);
                    Anew_double=double(Anews);
                    [Avec, Aeig]=eig(Anew_double);
                    Bnew_double=pinv(Avec,1e-14)*Bnews;
                    Cnew_double=Cnews*Avec;
                    for j=1:i
                        ws(j)=Bnew_double(j)*Cnew_double(j);
                    end
                    ps=diag(Aeig);
                    save(strjoin([pathname,'/s_', string(i)],''),'ps');
                    save(strjoin([pathname,'/w_', string(i)], ''),'ws');
                end
            end
            fd=fopen([pathname,'/readme.txt'],'w');
            t0='function=';
            t1='n=';t2='scale/n=';t3='digits=';
            t4='step=';t5='xmin=';t6='xmax=';
            t7='MRterms=';t8='Maxerror=';
            t9 = 'T='; t10 = 'Weight Function='; t11 = 'Weight Setting=';
            fprintf(fd,'%s',t0);fprintf(fd,'%s',handles.s);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t1);fprintf(fd,'%g',handles.n);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t9);fprintf(fd,'%g',handles.T);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t2);fprintf(fd,'%g',handles.scale);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t3);fprintf(fd,'%g',handles.digitss);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t4);fprintf(fd,'%g',handles.stepsize);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t5);fprintf(fd,'%g',handles.xmin);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t6);fprintf(fd,'%g',handles.xmax);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t10);fprintf(fd,'%s',handles.weight);fprintf(fd,'\r\n');
            fprintf(fd,'%s',t11);fprintf(fd,'%s',handles.weightset);fprintf(fd,'\r\n');
            fprintf(fd, '%s', t7);
            formatStr = [repmat('%g,', 1, numel(handles.mrterms)-1), '%g\r\n'];
            fprintf(fd, formatStr, handles.mrterms);
            fprintf(fd, '%s', t8);
            formatStr = [repmat('%g,', 1, numel(handles.errors)-1), '%g\r\n'];
            fprintf(fd, formatStr, handles.errors);
            fclose(fd);
            guidata(hObject, handles);
        end

        % Button pushed function: pushbutton6
        function pushbutton6_Callback(app, event)
            % --- Executes on button press in pushbutton6.
            
            % Create GUIDE-style callback args - Added by Migration Tool
            [hObject, eventdata, handles] = convertToGUIDECallbackArguments(app, event); %#ok<ASGLU>
            
            % hObject    handle to pushbutton6 (see GCBO)
            % eventdata  reserved - to be defined in a future version of MATLAB
            % handles    structure with handles and user data (see GUIDATA)
            handles.popindex=1;
            handles.err=[];
            handles.XX=0;
            handles.pp=[];
            handles.ww=[];
            handles.Anew=[];
            handles.Bnew=[];
            handles.Cnew=[];
            handles.radioindex=0;
            handles.radioindex2=0;
            handles.digitss=300;
            handles.xmax=0;
            handles.xmin=0;
            handles.n=0;
            handles.stepsize=0;
            handles.sclae=0;
            handles.errors=0;
            handles.errs=[];
            handles.s=[];
            handles.weight =[];
            handles.weightset =[];
            handles.T = 0;
            handles.accindex=1;
            handles.mrterms=0;
            handles.intindex=2;
            set(handles.pushbutton1,'Visible','Off');
            set(handles.pushbutton3,'Visible','Off');
            set(handles.pushbutton4,'Visible','Off');
            set(handles.pushbutton6,'Visible','Off');
            set(handles.text3,'Visible','Off');
            set(handles.text4,'Visible','Off');
            set(handles.text5,'Visible','Off');
            set(handles.edit1,'Visible','Off');
            set(handles.edit2,'Visible','Off');
            set(handles.edit3,'Visible','Off');
            set(handles.text7,'Visible','Off');
            set(handles.text8,'Visible','Off');
            set(handles.text11,'Visible','Off');
            set(handles.text13,'Visible','Off');
            set(handles.edit5,'Visible','Off');
            set(handles.text15,'Visible','Off');
            set(handles.edit6,'Visible','Off');
            set(handles.text16,'Visible','Off');
            set(handles.text19,'Visible','Off');
            set(handles.text20,'Visible','Off');
            set(handles.edit7,'Visible','Off');
            set(handles.edit8,'Visible','Off');
            set(handles.edit9,'Visible','Off');
            set(handles.edit10,'Visible','Off');
            set(handles.edit11,'Visible','Off');
            set(handles.text22,'Visible','Off');
            set(handles.text23,'Visible','Off');
            set(handles.text24,'Visible','Off');
            set(handles.text25,'Visible','Off');
            set(handles.text26,'Visible','Off');
            set(handles.popupmenu4,'Visible','Off');
            set(handles.popupmenu5,'Visible','Off');
            set(handles.popupmenu6,'Visible','Off');
            set(handles.text27,'Visible','Off');
            set(app.edit16,'Visible','Off');
            % Update handles structure
            guidata(hObject, handles);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create figure1 and hide until all components are created
            app.figure1 = uifigure('Visible', 'off');
            app.figure1.Position = [680 422 922 741];
            app.figure1.Name = 'VP-WBT';
            app.figure1.Resize = 'off';
            app.figure1.HandleVisibility = 'callback';
            app.figure1.Tag = 'figure1';

            % Create axes1
            app.axes1 = uiaxes(app.figure1);
            app.axes1.FontSize = 13.3333333333333;
            app.axes1.NextPlot = 'replace';
            app.axes1.Tag = 'axes1';
            app.axes1.Position = [425 195 430 404];

            % Create edit1
            app.edit1 = uieditfield(app.figure1, 'text');
            app.edit1.Tag = 'edit1';
            app.edit1.HorizontalAlignment = 'center';
            app.edit1.FontSize = 21.3333333333333;
            app.edit1.Position = [90 527 97 39];

            % Create text3
            app.text3 = uilabel(app.figure1);
            app.text3.Tag = 'text3';
            app.text3.HorizontalAlignment = 'center';
            app.text3.VerticalAlignment = 'top';
            app.text3.WordWrap = 'on';
            app.text3.FontSize = 26.6666666666667;
            app.text3.Position = [22 525 60 41];
            app.text3.Text = 'n';

            % Create text4
            app.text4 = uilabel(app.figure1);
            app.text4.Tag = 'text4';
            app.text4.HorizontalAlignment = 'center';
            app.text4.VerticalAlignment = 'top';
            app.text4.WordWrap = 'on';
            app.text4.FontSize = 26.6666666666667;
            app.text4.Position = [196 525 92 41];
            app.text4.Text = 'scale/n';

            % Create edit2
            app.edit2 = uieditfield(app.figure1, 'text');
            app.edit2.Tag = 'edit2';
            app.edit2.HorizontalAlignment = 'center';
            app.edit2.FontSize = 21.3333333333333;
            app.edit2.Position = [292 528 101 39];

            % Create edit3
            app.edit3 = uieditfield(app.figure1, 'text');
            app.edit3.Tag = 'edit3';
            app.edit3.HorizontalAlignment = 'center';
            app.edit3.FontSize = 21.3333333333333;
            app.edit3.Position = [292 436 102 39];

            % Create text5
            app.text5 = uilabel(app.figure1);
            app.text5.Tag = 'text5';
            app.text5.HorizontalAlignment = 'center';
            app.text5.VerticalAlignment = 'top';
            app.text5.WordWrap = 'on';
            app.text5.FontSize = 26.6666666666667;
            app.text5.Position = [205 432 66 41];
            app.text5.Text = 'max';

            % Create text6
            app.text6 = uilabel(app.figure1);
            app.text6.Tag = 'text6';
            app.text6.HorizontalAlignment = 'center';
            app.text6.VerticalAlignment = 'top';
            app.text6.WordWrap = 'on';
            app.text6.FontSize = 26.6666666666667;
            app.text6.Position = [8 671 117 33];
            app.text6.Text = 'Function';

            % Create edit4
            app.edit4 = uieditfield(app.figure1, 'text');
            app.edit4.Tag = 'edit4';
            app.edit4.HorizontalAlignment = 'center';
            app.edit4.FontSize = 21.3333333333333;
            app.edit4.Position = [133 666 487 39];

            % Create pushbutton1
            app.pushbutton1 = uibutton(app.figure1, 'push');
            app.pushbutton1.ButtonPushedFcn = createCallbackFcn(app, @pushbutton1_Callback, true);
            app.pushbutton1.Tag = 'pushbutton1';
            app.pushbutton1.FontSize = 21.3333333333333;
            app.pushbutton1.Position = [31 105 109 43];
            app.pushbutton1.Text = 'SOG';

            % Create text7
            app.text7 = uilabel(app.figure1);
            app.text7.Tag = 'text7';
            app.text7.HorizontalAlignment = 'center';
            app.text7.VerticalAlignment = 'top';
            app.text7.WordWrap = 'on';
            app.text7.FontSize = 26.6666666666667;
            app.text7.Position = [469 65 120 41];
            app.text7.Text = 'Maxerror';

            % Create text8
            app.text8 = uilabel(app.figure1);
            app.text8.Tag = 'text8';
            app.text8.VerticalAlignment = 'top';
            app.text8.WordWrap = 'on';
            app.text8.FontSize = 26.6666666666667;
            app.text8.Position = [606 65 198 41];
            app.text8.Text = '';

            % Create text9
            app.text9 = uilabel(app.figure1);
            app.text9.Tag = 'text9';
            app.text9.VerticalAlignment = 'top';
            app.text9.WordWrap = 'on';
            app.text9.FontSize = 10.6666666666667;
            app.text9.Position = [65 712 790 30];
            app.text9.Text = 'VPWBT-SOG(SOE)          Version 1.0 updated on 25.Feb.18th          Produced by Y.S.Lin Based on VPMR Prod. by J.Y.Liang and Z.X.Gao          GUI on MATLAB by Y.S.Lin and Z.X.Gao';

            % Create text10
            app.text10 = uilabel(app.figure1);
            app.text10.Tag = 'text10';
            app.text10.HorizontalAlignment = 'center';
            app.text10.VerticalAlignment = 'top';
            app.text10.WordWrap = 'on';
            app.text10.FontSize = 26.6666666666667;
            app.text10.Position = [20 567 94 41];
            app.text10.Text = 'Mode';

            % Create popupmenu1
            app.popupmenu1 = uidropdown(app.figure1);
            app.popupmenu1.Items = {'Select Mode', 'Pre-Experiment', 'Scale Experiment', 'Experiment', 'Test-Experiment', 'Pre-Experiment(Simple)', 'Scale Experiment(Simple)', 'Experiment(Simple)', 'Test-Experiment(Simple)'};
            app.popupmenu1.ValueChangedFcn = createCallbackFcn(app, @popupmenu1_Callback, true);
            app.popupmenu1.Tag = 'popupmenu1';
            app.popupmenu1.FontSize = 26.6666666666667;
            app.popupmenu1.BackgroundColor = [1 1 1];
            app.popupmenu1.Position = [133 572 254 40];
            app.popupmenu1.Value = 'Select Mode';

            % Create pushbutton3
            app.pushbutton3 = uibutton(app.figure1, 'push');
            app.pushbutton3.ButtonPushedFcn = createCallbackFcn(app, @pushbutton3_Callback, true);
            app.pushbutton3.Tag = 'pushbutton3';
            app.pushbutton3.FontSize = 21.3333333333333;
            app.pushbutton3.Position = [153 105 109 43];
            app.pushbutton3.Text = 'SOE';

            % Create pushbutton4
            app.pushbutton4 = uibutton(app.figure1, 'push');
            app.pushbutton4.ButtonPushedFcn = createCallbackFcn(app, @pushbutton4_Callback, true);
            app.pushbutton4.Tag = 'pushbutton4';
            app.pushbutton4.FontSize = 21.3333333333333;
            app.pushbutton4.Position = [59 18 293 43];
            app.pushbutton4.Text = 'Save results';

            % Create text11
            app.text11 = uilabel(app.figure1);
            app.text11.Tag = 'text11';
            app.text11.HorizontalAlignment = 'center';
            app.text11.VerticalAlignment = 'top';
            app.text11.WordWrap = 'on';
            app.text11.FontSize = 26.6666666666667;
            app.text11.Position = [600 121 204 33];
            app.text11.Text = ' ';

            % Create text12
            app.text12 = uilabel(app.figure1);
            app.text12.Tag = 'text12';
            app.text12.HorizontalAlignment = 'center';
            app.text12.VerticalAlignment = 'top';
            app.text12.WordWrap = 'on';
            app.text12.FontSize = 26.6666666666667;
            app.text12.Position = [468 157 376 33];
            app.text12.Text = ' ';

            % Create text13
            app.text13 = uilabel(app.figure1);
            app.text13.Tag = 'text13';
            app.text13.HorizontalAlignment = 'center';
            app.text13.VerticalAlignment = 'top';
            app.text13.WordWrap = 'on';
            app.text13.FontSize = 26.6666666666667;
            app.text13.Position = [25 382 108 41];
            app.text13.Text = 'MR terms';

            % Create edit5
            app.edit5 = uieditfield(app.figure1, 'text');
            app.edit5.Tag = 'edit5';
            app.edit5.HorizontalAlignment = 'center';
            app.edit5.FontSize = 21.3333333333333;
            app.edit5.Position = [143 386 251 39];

            % Create text15
            app.text15 = uilabel(app.figure1);
            app.text15.Tag = 'text15';
            app.text15.HorizontalAlignment = 'center';
            app.text15.VerticalAlignment = 'top';
            app.text15.WordWrap = 'on';
            app.text15.FontSize = 26.6666666666667;
            app.text15.Position = [356 26 128 33];
            app.text15.Text = 'Path';

            % Create edit6
            app.edit6 = uieditfield(app.figure1, 'text');
            app.edit6.Tag = 'edit6';
            app.edit6.HorizontalAlignment = 'center';
            app.edit6.FontSize = 21.3333333333333;
            app.edit6.Position = [483 20 291 39];

            % Create text16
            app.text16 = uilabel(app.figure1);
            app.text16.Tag = 'text16';
            app.text16.HorizontalAlignment = 'center';
            app.text16.VerticalAlignment = 'top';
            app.text16.WordWrap = 'on';
            app.text16.FontSize = 26.6666666666667;
            app.text16.Position = [22 478 60 41];
            app.text16.Text = 'digits';

            % Create edit7
            app.edit7 = uieditfield(app.figure1, 'text');
            app.edit7.Tag = 'edit7';
            app.edit7.HorizontalAlignment = 'center';
            app.edit7.FontSize = 21.3333333333333;
            app.edit7.Position = [89 482 98 39];

            % Create edit8
            app.edit8 = uieditfield(app.figure1, 'text');
            app.edit8.Tag = 'edit8';
            app.edit8.HorizontalAlignment = 'center';
            app.edit8.FontSize = 21.3333333333333;
            app.edit8.Position = [89 434 97 39];

            % Create text19
            app.text19 = uilabel(app.figure1);
            app.text19.Tag = 'text19';
            app.text19.HorizontalAlignment = 'center';
            app.text19.VerticalAlignment = 'top';
            app.text19.WordWrap = 'on';
            app.text19.FontSize = 26.6666666666667;
            app.text19.Position = [15 432 69 41];
            app.text19.Text = 'min';

            % Create edit9
            app.edit9 = uieditfield(app.figure1, 'text');
            app.edit9.Tag = 'edit9';
            app.edit9.HorizontalAlignment = 'center';
            app.edit9.FontSize = 21.3333333333333;
            app.edit9.Position = [292 483 101 39];

            % Create text20
            app.text20 = uilabel(app.figure1);
            app.text20.Tag = 'text20';
            app.text20.HorizontalAlignment = 'center';
            app.text20.VerticalAlignment = 'top';
            app.text20.WordWrap = 'on';
            app.text20.FontSize = 26.6666666666667;
            app.text20.Position = [198 480 81 41];
            app.text20.Text = 'step';

            % Create edit10
            app.edit10 = uieditfield(app.figure1, 'text');
            app.edit10.Tag = 'edit10';
            app.edit10.HorizontalAlignment = 'center';
            app.edit10.FontSize = 21.3333333333333;
            app.edit10.Position = [262 619 316 39];

            % Create text22
            app.text22 = uilabel(app.figure1);
            app.text22.Tag = 'text22';
            app.text22.HorizontalAlignment = 'center';
            app.text22.VerticalAlignment = 'top';
            app.text22.WordWrap = 'on';
            app.text22.FontSize = 26.6666666666667;
            app.text22.Position = [572 618 94 41];
            app.text22.Text = 'SP';

            % Create edit11
            app.edit11 = uieditfield(app.figure1, 'text');
            app.edit11.Tag = 'edit11';
            app.edit11.HorizontalAlignment = 'center';
            app.edit11.FontSize = 21.3333333333333;
            app.edit11.Position = [665 620 178 39];

            % Create text23
            app.text23 = uilabel(app.figure1);
            app.text23.Tag = 'text23';
            app.text23.HorizontalAlignment = 'center';
            app.text23.VerticalAlignment = 'top';
            app.text23.WordWrap = 'on';
            app.text23.FontSize = 26.6666666666667;
            app.text23.Position = [107 478 60 41];
            app.text23.Text = '300';

            % Create text24
            app.text24 = uilabel(app.figure1);
            app.text24.Tag = 'text24';
            app.text24.HorizontalAlignment = 'center';
            app.text24.VerticalAlignment = 'top';
            app.text24.WordWrap = 'on';
            app.text24.FontSize = 26.6666666666667;
            app.text24.Position = [295 478 95 41];
            app.text24.Text = 'Default';

            % Create text25
            app.text25 = uilabel(app.figure1);
            app.text25.Tag = 'text25';
            app.text25.HorizontalAlignment = 'center';
            app.text25.VerticalAlignment = 'top';
            app.text25.WordWrap = 'on';
            app.text25.FontSize = 26.6666666666667;
            app.text25.Position = [107 432 60 41];
            app.text25.Text = '0';

            % Create text26
            app.text26 = uilabel(app.figure1);
            app.text26.Tag = 'text26';
            app.text26.HorizontalAlignment = 'center';
            app.text26.VerticalAlignment = 'top';
            app.text26.WordWrap = 'on';
            app.text26.FontSize = 26.6666666666667;
            app.text26.Position = [25 304 94 33];
            app.text26.Text = 'Integral';

            % Create popupmenu4
            app.popupmenu4 = uidropdown(app.figure1);
            app.popupmenu4.Items = {'Select Method', 'Quadgk', 'Gauss base points'};
            app.popupmenu4.ValueChangedFcn = createCallbackFcn(app, @popupmenu4_Callback, true);
            app.popupmenu4.Tag = 'popupmenu4';
            app.popupmenu4.FontSize = 26.6666666666667;
            app.popupmenu4.BackgroundColor = [1 1 1];
            app.popupmenu4.Position = [133 300 254 40];
            app.popupmenu4.Value = 'Select Method';

            % Create text27
            app.text27 = uilabel(app.figure1);
            app.text27.Tag = 'text27';
            app.text27.HorizontalAlignment = 'center';
            app.text27.VerticalAlignment = 'top';
            app.text27.WordWrap = 'on';
            app.text27.FontSize = 26.6666666666667;
            app.text27.Position = [24 258 100 33];
            app.text27.Text = 'Accuracy';

            % Create popupmenu5
            app.popupmenu5 = uidropdown(app.figure1);
            app.popupmenu5.Items = {'Select Integral Accuracy', 'Low(1e-8)', 'Mid(1e-10)', 'High(1e-12)', 'High+(1e-14)'};
            app.popupmenu5.ValueChangedFcn = createCallbackFcn(app, @popupmenu5_Callback, true);
            app.popupmenu5.Tag = 'popupmenu5';
            app.popupmenu5.FontSize = 26.6666666666667;
            app.popupmenu5.BackgroundColor = [1 1 1];
            app.popupmenu5.Position = [133 255 255 40];
            app.popupmenu5.Value = 'Select Integral Accuracy';

            % Create popupmenu6
            app.popupmenu6 = uidropdown(app.figure1);
            app.popupmenu6.Items = {'Select Integral Accuracy', 'Low(5e6)', 'Mid(1e7)', 'High(5e7)', 'High+(1e8)'};
            app.popupmenu6.ValueChangedFcn = createCallbackFcn(app, @popupmenu6_Callback, true);
            app.popupmenu6.Tag = 'popupmenu6';
            app.popupmenu6.FontSize = 26.6666666666667;
            app.popupmenu6.BackgroundColor = [1 1 1];
            app.popupmenu6.Position = [133 255 255 40];
            app.popupmenu6.Value = 'Select Integral Accuracy';

            % Create pushbutton6
            app.pushbutton6 = uibutton(app.figure1, 'push');
            app.pushbutton6.ButtonPushedFcn = createCallbackFcn(app, @pushbutton6_Callback, true);
            app.pushbutton6.Tag = 'pushbutton6';
            app.pushbutton6.FontSize = 21.3333333333333;
            app.pushbutton6.Position = [274 105 110 43];
            app.pushbutton6.Text = 'Reset';

            % Create text28
            app.text28 = uilabel(app.figure1);
            app.text28.Tag = 'text28';
            app.text28.HorizontalAlignment = 'center';
            app.text28.VerticalAlignment = 'top';
            app.text28.WordWrap = 'on';
            app.text28.FontSize = 26.6666666666667;
            app.text28.Position = [0 383 89 41];
            app.text28.Text = 'sc-min';

            % Create text29
            app.text29 = uilabel(app.figure1);
            app.text29.Tag = 'text29';
            app.text29.HorizontalAlignment = 'center';
            app.text29.VerticalAlignment = 'top';
            app.text29.WordWrap = 'on';
            app.text29.FontSize = 26.6666666666667;
            app.text29.Position = [194 390 94 35];
            app.text29.Text = 'sc-max';

            % Create edit13
            app.edit13 = uieditfield(app.figure1, 'text');
            app.edit13.Tag = 'edit13';
            app.edit13.HorizontalAlignment = 'center';
            app.edit13.FontSize = 21.3333333333333;
            app.edit13.Position = [89 386 99 39];

            % Create edit14
            app.edit14 = uieditfield(app.figure1, 'text');
            app.edit14.Tag = 'edit14';
            app.edit14.HorizontalAlignment = 'center';
            app.edit14.FontSize = 21.3333333333333;
            app.edit14.Position = [292 386 101 39];

            % Create text30
            app.text30 = uilabel(app.figure1);
            app.text30.Tag = 'text30';
            app.text30.HorizontalAlignment = 'center';
            app.text30.VerticalAlignment = 'top';
            app.text30.WordWrap = 'on';
            app.text30.FontSize = 26.6666666666667;
            app.text30.Position = [0 338 90 41];
            app.text30.Text = 'sc-step';

            % Create edit15
            app.edit15 = uieditfield(app.figure1, 'text');
            app.edit15.Tag = 'edit15';
            app.edit15.HorizontalAlignment = 'center';
            app.edit15.FontSize = 21.3333333333333;
            app.edit15.Position = [89 343 98 39];

            % Create popupmenu7
            app.popupmenu7 = uidropdown(app.figure1);
            app.popupmenu7.Items = {'Record'};
            app.popupmenu7.ValueChangedFcn = createCallbackFcn(app, @popupmenu7_Callback, true);
            app.popupmenu7.Tag = 'popupmenu7';
            app.popupmenu7.FontSize = 26.6666666666667;
            app.popupmenu7.BackgroundColor = [1 1 1];
            app.popupmenu7.Position = [620 665 223 40];
            app.popupmenu7.Value = 'Record';

            % Create TTextAreaLabel
            app.TTextAreaLabel = uilabel(app.figure1);
            app.TTextAreaLabel.HorizontalAlignment = 'right';
            app.TTextAreaLabel.FontSize = 24;
            app.TTextAreaLabel.Position = [8 626 25 32];
            app.TTextAreaLabel.Text = 'T';

            % Create edit16
            app.edit16 = uitextarea(app.figure1);
            app.edit16.FontSize = 24;
            app.edit16.Position = [48 619 70 41];
            app.edit16.Value = {'1024'};

            % Create WeightDropDownLabel
            app.WeightDropDownLabel = uilabel(app.figure1);
            app.WeightDropDownLabel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.WeightDropDownLabel.HorizontalAlignment = 'center';
            app.WeightDropDownLabel.FontSize = 26;
            app.WeightDropDownLabel.Position = [26 209 85 34];
            app.WeightDropDownLabel.Text = 'Weight';

            % Create WeightDropDown
            app.WeightDropDown = uidropdown(app.figure1);
            app.WeightDropDown.Items = {'Select Weight', '1', '(x+a).^(-1/2)', 'exp(-s.*x)', 'Custom'};
            app.WeightDropDown.FontSize = 26;
            app.WeightDropDown.BackgroundColor = [1 1 1];
            app.WeightDropDown.Position = [133 209 254 34];
            app.WeightDropDown.Value = 'Select Weight';

            % Create WeightSettingEditFieldLabel
            app.WeightSettingEditFieldLabel = uilabel(app.figure1);
            app.WeightSettingEditFieldLabel.HorizontalAlignment = 'right';
            app.WeightSettingEditFieldLabel.FontSize = 24;
            app.WeightSettingEditFieldLabel.Position = [25 161 161 32];
            app.WeightSettingEditFieldLabel.Text = 'Weight Setting';

            % Create WeightSettingEditField
            app.WeightSettingEditField = uieditfield(app.figure1, 'text');
            app.WeightSettingEditField.FontSize = 24;
            app.WeightSettingEditField.Position = [198 158 188 38];

            % Show the figure after all components are created
            app.figure1.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = VP_WBT

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.figure1)

                % Execute the startup function
                runStartupFcn(app, @VP_WBT_OpeningFcn)
            else

                % Focus the running singleton app
                figure(runningApp.figure1)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.figure1)
        end
    end
end