classdef parse_LF_data<dynamicprops
    properties
        frames=[];
        calib_eV=[];
        calib_nm=[];
        file_spe=[];
        file_photons=[];
    end
    
    methods
        
        %% CONSTRUCTOR METHOD
        function obj=parse_LF_data(file_path_spe,file_path_photon_stream,ht3_file_lifetime,buffer_size,pico_path)
            
            % Creates a SPEreader object that is wrapped into the analyis object.
            % then reads the reader object and sabed
            
            if nargin<2
                file_path_photon_stream=0;
            end
            
            obj.file_spe=file_path_spe; % store as property
            obj.file_photons=file_path_photon_stream;
            reader_object=SpeReader(file_path_spe); % create SpeReader object
            
            if ~isprop(obj,'SPE_reader') % add property to the class
                dummy=addprop(obj,'SPE_reader');
            end
            obj.SPE_reader=reader_object; % write the SpeReader object to the class as property
            
            if ~isprop(obj,'frames')% add property to the class
                dummy=addprop(obj,'frames');
            end
            
            if ~isprop(obj,'frameTimes')% add property to the class
                dummy=addprop(obj,'frameTimes');
            end
            
            [frameTimes,frames]=read(reader_object); % save the frame array as created by reading the reader object to class as property.
            obj.frames=frames; % is the spectral matrix Pixels X numFrames.
            obj.frameTimes=frameTimes; % is the spectral matrix Pixels X numFrames.
            
            if file_path_photon_stream~=0
                photons_object=photons(file_path_photon_stream,buffer_size,pico_path)
                
                if ~isprop(obj,'photon_stream')% add property to the class
                    dummy=addprop(obj,'photon_stream');
                end
                obj.photon_stream=photons_object;
                
                photons_lifetime=photons(ht3_file_lifetime,buffer_size,pico_path)
                if ~isprop(obj,'photons_lifetime')% add property to the class
                    dummy=addprop(obj,'photons_lifetime');
                end
                obj.photons_lifetime=photons_lifetime;
            end
            % obtain calibration file and perform simple analysis
            obj.get_calib()
            obj.simple_analysis()
            obj.plot_simple_analysis()
        end
        
        %% other methods to parse data.
        function obj=get_calib(obj)
            %extract calibration from the SPE XLM.
            %Store the XML data in a temp *.xml file
            filename = ['temp' '.xml'];
            fid = fopen(filename,'Wt');
            fwrite(fid,obj.SPE_reader.FooterXML);
            fclose(fid);
            
            % Read the file into an XML model object
            xmlstruct = xml2struct(filename);
            
            % Delete the temp file
            delete(filename);
            
            % Parse the XML model object
            if ~isprop(obj,'footerxmlstruct')% add property to the class
                dummy=addprop(obj,'footerxmlstruct');
            end
            obj.footerxmlstruct=xmlstruct;
            
            %%get_calibration
            fid = fopen('temp_calib.txt','wt');
            fprintf(fid,obj.footerxmlstruct.SpeFormat.Calibrations.WavelengthMapping.Wavelength.Text);
            fclose(fid);
            
            calib=csvread('temp_calib.txt')
            delete('temp_calib.txt')
            
            if ~isprop(obj,'calib_nm')% add property to the class
                dummy=addprop(obj,'calib_nm');
            end
            if ~isprop(obj,'calib_eV')% add property to the class
                dummy=addprop(obj,'calib_eV');
            end
            
            obj.calib_nm=calib;
            obj.calib_eV=1234./calib*1000;
        end
        function obj=get_intensity_trace_HH(obj,bin_width_ps)
            obj.photon_stream.picoquant_bin()
            obj.photon_stream.get_intensity_w_SPEmarkers(obj.file_photons(1:end-4),'intensity',bin_width_ps)
        end
        
        
        %% time series exploration and analysis.
        function obj=truncate_frames_and_calib(obj,lower_pixel,upper_pixel)
            % saves part of the the frame array as a new frame array to get rid of
            % laser leakthrough etc.
            obj.frames=obj.frames(lower_pixel:upper_pixel,:);
            obj.calib_eV=obj.calib_eV(lower_pixel:upper_pixel);
            
        end
        function obj=baseline_correct(obj,baseline)
            for i=1:length(obj.frames(1,:))
                obj.frames(:,i)=obj.frames(:,i)-baseline;
            end
        end
        function obj=plot_frame(obj,lower_frame,upper_frame)
            if nargin<3
                figure()
                plot(obj.calib_eV,obj.frames(:,lower_frame))
                title(strcat('Frame','  ', num2str(lower_frame)))
                
            else
                figure()
                plot(obj.calib_eV,mean(obj.frames(:,lower_frame:upper_frame),2))
                title(strcat('Frames','  ', num2str(lower_frame),' - ',num2str(upper_frame)))
                
                
            end
            
        end
        function obj=simple_analysis(obj,plot)
            
            iterations=obj.SPE_reader.NumberOfFrames
            
            intensity=zeros(1,iterations);
            sum_spectrum_sqaured=zeros(1,iterations);
            max_energy=zeros(1,iterations);
            max_intensity=zeros(1,iterations);
            
            for i=1:iterations;
                intensity(i)=sum(obj.frames(:,i));
                sum_spectrum_sqaured(i)=sum(obj.frames(:,i).^2)/intensity(i).^2;
                [max_intensity(i),max_energy(i)]=max(obj.frames(:,i));
            end
            
            % do some histogramming of the time trace properties.
            
            
            
            
            if ~isprop(obj,'simple_analysis_results')% add property to the class
                dummy=addprop(obj,'simple_analysis_results');
            end
            
            obj.simple_analysis_results=struct('intensity',intensity,'sum_spectrum_sqaured',sum_spectrum_sqaured,'max_energy',max_energy,'max_intensity',max_intensity); % save the frame array as created by reading the reader object to class as property.
        end
        function obj=plot_simple_analysis(obj)
            figure('position', [0, 0, 1400, 900])
            title(obj.file_spe)
            subplot(5,1,1)
            plot(obj.simple_analysis_results.intensity)
            title('Intensity Camera')
            subplot(5,1,2)
            plot(obj.simple_analysis_results.sum_spectrum_sqaured./obj.simple_analysis_results.intensity)
            title('(Sum Spectrum^2)/Intensity^2 Camera')
            subplot(5,1,3)
            plot(obj.calib_eV(obj.simple_analysis_results.max_energy))
            title('Max Energy')
            subplot(5,1,4)
            plot(obj.simple_analysis_results.max_intensity)
            title('Max Intensity Camera')
            subplot(5,1,5)
            contourf(obj.frames, 'Linecolor','none')
        end
        function obj=contour_spectra_photons(obj,ylim_spectra)
            %% plot contour of spectra together with intensity trace for the photon-stream of HH.
            t_0=obj.photon_stream.intensity.marker_times(1);
            
            %get histogram of the photon counts and average spectrum.
            intensity=sum(obj.photon_stream.intensity.trace,2)./(obj.photon_stream.intensity.bin_width/1E12);
            h=histogram(intensity,50);
            average_spec=mean(obj.frames,2)
            
            figure('pos',[200 200 1000 500])
            subplot(2,2,1)
            
            contourf(obj.photon_stream.intensity.marker_times(1:end-1)-t_0,obj.calib_eV,obj.frames, 'Linecolor','none')
            colormap('jet')
            xlim([0,max(obj.photon_stream.intensity.marker_times)-t_0])
            ylabel('Energy [meV]')
            xlabel('time [s]')
            ylim(ylim_spectra)
            title(obj.photon_stream.fname,'Interpreter', 'none')
            set(gca,'fontsize',14)
            
            subplot(2,2,2)
            plot(average_spec,obj.calib_eV,'linewidth',1,'color','blue')
            ylim(ylim_spectra)
            xlim([0,max(average_spec)+max(average_spec)/4])
            set(gca,'fontsize',14)
            xlabel('Intensity')
            ylabel('Energy [meV]')
            
            
            subplot(2,2,3)
            plot(obj.photon_stream.time-t_0,sum(obj.photon_stream.intensity.trace,2)./(obj.photon_stream.intensity.bin_width/1E12),'Color','blue')
            hold on
            plot(obj.photon_stream.intensity.marker_times(2:end)-t_0,obj.photon_stream.intensity.marker_trace(2:end)./2,'Linewidth',2,'Color','red')
            xlim([0,max(obj.photon_stream.intensity.marker_times)-t_0])
            foo=max(sum(obj.photon_stream.intensity.trace,2)./(obj.photon_stream.intensity.bin_width/1E12));
            ylim([0,(foo+foo/4)])
            ylabel('cps')
            xlabel('time [s]')
            set(gca,'fontsize',14)

            subplot(2,2,4)
            plot(h.Values,h.BinEdges(1:end-1),'linewidth',1,'color','blue')
            ylim([0,(foo+foo/4)])
            set(gca,'fontsize',14)
            ylabel('cps')
            xlabel('occurences')

            saveas(gcf,strcat('contour_photons',obj.photon_stream.fname(1:end-4)),'svg')

        end
        function obj=contour_spectra(obj,ylimits)
            %% plot contour of spectra without any photon-stream data.
            figure('pos',[100 100 1000 300])
            contourf([1:1:size(obj.frames,2)],obj.calib_eV,obj.frames, 'Linecolor','none')
            colormap('jet')
            ylabel('Energy [meV]')
            xlabel('Frame')
            ylim(ylimits)
            set(gca,'FontSize', 15)
            title(obj.file_spe,'interpreter','none')
            fname=char(obj.file_spe);
            set(gca,'fontsize',17)
            saveas(gcf,strcat('contour',fname(1:end-4)),'svg')
        end
        
        function obj=correlations(obj)
            %% plot correlation maps
            max_energy=obj.calib_eV(obj.simple_analysis_results.max_energy(1:end-1));
            intensity=obj.photon_stream.intensity.marker_trace(2:end);
            
            figure()
            title(obj.photon_stream.fname)
            subplot(2,2,1)
            plot(intensity,max_energy,'o');
            ylabel('Emision Maximum [eV]')
            xlabel('Intensity [cps]')
            subplot(2,2,2)
            plot(obj.simple_analysis_results.max_intensity,obj.calib_eV(obj.simple_analysis_results.max_energy),'o','Color','green')
            ylabel('Emision Maximum [eV]')
            xlabel('Maximum Camera Intensity [cps]')
            subplot(2,2,3)
            plot(obj.simple_analysis_results.sum_spectrum_sqaured(1:end-1)./(obj.photon_stream.intensity.marker_trace(2:end)).^2,obj.calib_eV(obj.simple_analysis_results.max_energy(1:end-1)),'o','Color','red')
            ylabel('Emission Maximum [eV]')
            xlabel('(Sum Spectrum)^2/Intensity^2')
            
        end
        function obj=create_spectral_series_video(obj,vid_fname)
            % vid_fname needs to be string with .avi ending.
            v = VideoWriter(vid_fname);
            open(v);
            
            figure()
            plot(obj.frames(:,1))
            
            axis tight
            set(gca,'nextplot','replacechildren');
            
            for k = 1:length(obj.frames(1,:))
                plot(obj.calib_eV,obj.frames(:,k))
                ylim=get(gca,'ylim');
                xlim=get(gca,'xlim')
                text(xlim(1)+xlim(1)/50,ylim(2)-ylim(2)/10,strcat('Frame ',' ',num2str(k)),'FontSize',14)
                frame = getframe(gcf);
                writeVideo(v,frame);
            end
            
            close(v);
        end
        function obj=presort_frames(obj,params_for_model_spectrum);
            %% compares the spectral time-series to a frame created by the model to dismiss really bad frames.
            %% get a reasonable frame with the model.
            
            lscerr=[];
            
            model_frame=obj.lineshape_model(obj.calib_eV,params_for_model_spectrum);
            [model_energy,model_norm]=obj.center_spectrum(obj.calib_eV,model_frame);
            
            %interpolate the centered model frame
            model_cent_interp=interp1(model_energy,model_norm,linspace(-200,200,500));
            
            figure()
            plot(linspace(-200,200,500),model_cent_interp)
            title('model spectrum used for presorting')
            xlabel('Energy [meV]')
            ylabel('Spectrum')
            
            for k=1:length(obj.frames(1,:));
                
                [cent_energy,cent_spec]=obj.center_spectrum(obj.calib_eV,obj.frames(:,k));
                spectrum_interp=interp1(cent_energy,cent_spec,linspace(-200,200,500));
                residual=(model_cent_interp-spectrum_interp).^2;
                residual(isnan(residual))=0;
                lscerr(k)=sum(residual);
                
            end
            
            if ~isprop(obj,'presort_error')% add property to the class
                dummy=addprop(obj,'presort_error');
            end
            
            obj.presort_error=lscerr;
            
            figure()
            plot(lscerr)
            xlabel('frame')
            ylabel('residual model frame')
            
        end
        
        %% create spectrum property.
        function obj=create_spectrum_struct(obj,energy,spectrum,spectrum_id);
            
            if ~isprop(obj,spectrum_id)% add property to the class
                dummy=addprop(obj,spectrum_id);
            end
            
            obj.(spectrum_id)=struct('energy',energy,'spectrum',spectrum);
            
        end
        function obj=create_spectrum_struct_average(obj,framevector,spectrum_id);
            
            energy=obj.calib_eV;
            spectrum=mean(obj.frames(:,framevector),2);
            
            if ~isprop(obj,spectrum_id)% add property to the class
                dummy=addprop(obj,spectrum_id);
            end
            
            obj.(spectrum_id)=struct('energy',energy,'spectrum',spectrum);
            
        end
        function obj=create_spectra_structs_from_presorting(obj,presort_threshold);
            framevector=[];
            for k=1:length(obj.frames(1,:));
                spectrum_id=strcat('frame',num2str(k));
                if obj.presort_error(k)<presort_threshold
                    obj.create_spectrum_struct(obj.calib_eV,obj.frames(:,k),spectrum_id)
                    framevector=[framevector,k];
                end
            end
            if ~isprop(obj,'framevector')% add property to the class
                dummy=addprop(obj,'framevector');
                
            end
            obj.framevector=framevector;
        end
 
        %% different fitting procedures.
        function obj=fit_lineshape_prefac_E0_HWHM(obj,spectrum_id,parameters,n_multistarts,fit_name);
            
            %retrive spectrum from the spectrum_id
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            
            %center and normalize the spectrum
            [center_energy,center_spectrum]=center_spec(energy,spectrum);
            
            %% set parameters for fit
            prefactor=parameters(1);
            E0=parameters(2);
            E_fine=parameters(3);
            T=parameters(4);
            E1=parameters(5);
            E2=parameters(6);
            E3=parameters(7);
            S1=parameters(8);
            S2=parameters(9);
            S3=parameters(10);
            HWHM_0=parameters(11);
            HWHM_T=parameters(12);
            
            fit_func=@(params,cent_energy)obj.lineshape_model(center_energy,[params(2),params(3),E_fine,T,E1,E2,E3,S1,S2,S3,params(1),HWHM_T]);
            
            %set lower and upper bound
            lb = [0,0,-200];
            ub = [100,1.1,200];
            
            %set arbitrary intital guess
            p0=[10,1,0];
            %  opts=optimoptions(@lsqcurvefit,'StepTolerance',1E-10,'FunctionTolerance',1E-10);
            
            problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',fit_func,'lb',lb,'ub',ub,'xdata',double(center_energy),'ydata',double(center_spectrum));
            ms = MultiStart('PlotFcns',@gsplotbestf);
            [fit_output,error_out] = run(ms,problem,n_multistarts)
            
            fit_params=[fit_output(2),fit_output(3),E_fine,T,E1,E2,E3,S1,S2,S3,fit_output(1),HWHM_T];
            errormulti=error_out;
            
            %obtain the FWHM for the underlying spectrum.
            FWHM=obj.get_FWHM_model(center_energy,fit_params);
            
            % adding the model parametes to the struct of the specific
            % spectrum.
            [obj.(spectrum_id).(fit_name).fit_params]=fit_params;
            [obj.(spectrum_id).(fit_name).fit_error]=errormulti;
            [obj.(spectrum_id).(fit_name).FWHM]=FWHM;
        end
        function obj=fit_lineshape_prefac_E0_HWHM_wo_center(obj,spectrum_id,parameters,n_multistarts,fit_name);
            
            %retrive spectrum from the spectrum_id
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            
            %center and normalize the spectrum
            %[center_energy,center_spectrum]=center_spec(energy,spectrum);
            center_energy=energy;
            center_spectrum=spectrum/max(spectrum);
            %% set parameters for fit
            prefactor=parameters(1);
            E0=parameters(2);
            E_fine=parameters(3);
            T=parameters(4);
            E1=parameters(5);
            E2=parameters(6);
            E3=parameters(7);
            S1=parameters(8);
            S2=parameters(9);
            S3=parameters(10);
            HWHM_0=parameters(11);
            HWHM_T=parameters(12);
            
            fit_func=@(params,cent_energy)obj.lineshape_model(center_energy,[params(2),params(3),E_fine,T,E1,E2,E3,S1,S2,S3,params(1),HWHM_T]);
            
            %set lower and upper bound
            lb = [0,0,-3000];
            ub = [300,1.1,2040];
            
            %set arbitrary intital guess
            p0=[10,1,2000];
            %  opts=optimoptions(@lsqcurvefit,'StepTolerance',1E-10,'FunctionTolerance',1E-10);
            
            problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',fit_func,'lb',lb,'ub',ub,'xdata',double(center_energy),'ydata',double(center_spectrum));
            ms = MultiStart('PlotFcns',@gsplotbestf);
            [fit_output,error_out] = run(ms,problem,n_multistarts)
            
            fit_params=[fit_output(2),fit_output(3),E_fine,T,E1,E2,E3,S1,S2,S3,fit_output(1),HWHM_T];
            errormulti=error_out;
            
            %obtain the FWHM for the underlying spectrum.
            FWHM=obj.get_FWHM_model(center_energy,fit_params);
            
            % adding the model parametes to the struct of the specific
            % spectrum.
            [obj.(spectrum_id).(fit_name).fit_params]=fit_params;
            [obj.(spectrum_id).(fit_name).fit_error]=errormulti;
            [obj.(spectrum_id).(fit_name).FWHM]=FWHM;
        end
        function obj=fit_lineshape_prefac_HWHM(obj,spectrum_id,parameters,n_multistarts,fit_name);
            
            %retrive spectrum from the spectrum_id
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            
            %center and normalize the spectrum
            [center_energy,center_spectrum]=center_spec(energy,spectrum);
            
            %% set parameters for fit
            prefactor=parameters(1);
            E0=parameters(2);
            E_fine=parameters(3);
            T=parameters(4);
            E1=parameters(5);
            E2=parameters(6);
            E3=parameters(7);
            S1=parameters(8);
            S2=parameters(9);
            S3=parameters(10);
            HWHM_0=parameters(11);
            HWHM_T=parameters(12);
            
            fit_func=@(params,cent_energy)obj.lineshape_model(center_energy,[params(2),E0,E_fine,T,E1,E2,E3,S1,S2,S3,params(1),HWHM_T]);
            
            %set lower and upper bound
            lb = [0,0];
            ub = [100,1.1];
            
            %set arbitrary intital guess
            p0=[10,1];
            %  opts=optimoptions(@lsqcurvefit,'StepTolerance',1E-10,'FunctionTolerance',1E-10);
            
            problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',fit_func,'lb',lb,'ub',ub,'xdata',double(center_energy),'ydata',double(center_spectrum));
            ms = MultiStart('PlotFcns',@gsplotbestf);
            [fit_output,error_out] = run(ms,problem,n_multistarts)
            
            fit_params=[fit_output(2),E0,E_fine,T,E1,E2,E3,S1,S2,S3,fit_output(1),HWHM_T];
            errormulti=error_out;
            
            %obtain the FWHM for the underlying spectrum.
            FWHM=obj.get_FWHM_model(center_energy,fit_params);
            
            % adding the model parametes to the struct of the specific
            % spectrum.
            [obj.(spectrum_id).(fit_name).fit_params]=fit_params;
            [obj.(spectrum_id).(fit_name).fit_error]=errormulti;
            [obj.(spectrum_id).(fit_name).FWHM]=FWHM;
        end
        function obj=fit_lineshape_HWHM(obj,spectrum_id,parameters,n_multistarts,fit_name);
            
            %retrive spectrum from the spectrum_id
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            
            %center and normalize the spectrum
            [center_energy,center_spectrum]=center_spec(energy,spectrum);
            
            %% set parameters for fit
            prefactor=parameters(1);
            E0=parameters(2);
            E_fine=parameters(3);
            T=parameters(4);
            E1=parameters(5);
            E2=parameters(6);
            E3=parameters(7);
            S1=parameters(8);
            S2=parameters(9);
            S3=parameters(10);
            HWHM_0=parameters(11);
            HWHM_T=parameters(12);
            
            fit_func=@(params,cent_energy)obj.lineshape_model(center_energy,[prefactor,E0,E_fine,T,E1,E2,E3,S1,S2,S3,params(1),HWHM_T]);
            
            %set lower and upper bound
            lb = [0];
            ub = [100];
            
            %set arbitrary intital guess
            p0=[10];
            %  opts=optimoptions(@lsqcurvefit,'StepTolerance',1E-10,'FunctionTolerance',1E-10);
            
            problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',fit_func,'lb',lb,'ub',ub,'xdata',double(center_energy),'ydata',double(center_spectrum));
            ms = MultiStart('PlotFcns',@gsplotbestf);
            [fit_output,error_out] = run(ms,problem,n_multistarts)
            
            fit_params=[prefactor,E0,E_fine,T,E1,E2,E3,S1,S2,S3,fit_output(1),HWHM_T];
            errormulti=error_out;
            
            %obtain the FWHM for the underlying spectrum.
            FWHM=obj.get_FWHM_model(center_energy,fit_params);
            
            % adding the model parametes to the struct of the specific
            % spectrum.
            [obj.(spectrum_id).(fit_name).fit_params]=fit_params;
            [obj.(spectrum_id).(fit_name).fit_error]=errormulti;
            [obj.(spectrum_id).(fit_name).FWHM]=FWHM;
        end
        function obj=fit_lineshape_prefac_E0_S1_S2_S3_E3_HWHM(obj,spectrum_id,parameters,n_multistarts,fit_name);
            
            %retrive spectrum from the spectrum_id
            disp(spectrum_id)
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            
            %center and normalize the spectrum
            [center_energy,center_spectrum]=center_spec(energy,spectrum);
            
            %% set parameters for fit
            prefactor=parameters(1);
            E0=parameters(2);
            E_fine=parameters(3);
            T=parameters(4);
            E1=parameters(5);
            E2=parameters(6);
            E3=parameters(7);
            S1=parameters(8);
            S2=parameters(9);
            S3=parameters(10);
            HWHM_0=parameters(11);
            HWHM_T=parameters(12);
            
            fit_func=@(params,cent_energy)obj.lineshape_model(cent_energy,[params(7),params(6),E_fine,T,E1,E2,params(4),params(1),params(2),params(5),params(3),HWHM_T]);
            
            % % %set lower and upper bound
            lb = [0,0,0,0,0,-100,0];
            ub = [10,10,100,10,10,100,1.5];
            %
            % %set arbitrary intital guess
            p0=[0.2,0.2,2,4,0.9,0,1];
            %            opts=optimoptions(@lsqcurvefit,'StepTolerance',1E-10,'FunctionTolerance',1E-10);
            
            problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',fit_func,'lb',lb,'ub',ub,'xdata',double(center_energy),'ydata',double(center_spectrum));
            ms = MultiStart('PlotFcns',@gsplotbestf);
            [fit_output,error_out] = run(ms,problem,n_multistarts)
            
            fit_params=[fit_output(7),fit_output(6),E_fine,T,E1,E2,fit_output(4),fit_output(1),fit_output(2),fit_output(5),fit_output(3),HWHM_T];
            errormulti=error_out;
            
            %obtain the FWHM for the underlying spectrum.
            FWHM=obj.get_FWHM_model(center_energy,fit_params);
            
            % adding the model parametes to the struct of the specific
            % spectrum.
            [obj.(spectrum_id).(fit_name).fit_params]=fit_params;
            [obj.(spectrum_id).(fit_name).fit_error]=errormulti;
            [obj.(spectrum_id).(fit_name).FWHM]=FWHM;
        end
        function obj=fit_lineshape_prefac_E0_S1_S2_S3_E3_HWHM_wo_center(obj,spectrum_id,parameters,n_multistarts,fit_name);
            
            %retrive spectrum from the spectrum_id
            disp(spectrum_id)
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            
            %center and normalize the spectrum
            %  [center_energy,center_spectrum]=center_spec(energy,spectrum);
            center_energy=energy;
            center_spectrum=spectrum/max(spectrum);
            %% set parameters for fit
            prefactor=parameters(1);
            E0=parameters(2);
            E_fine=parameters(3);
            T=parameters(4);
            E1=parameters(5);
            E2=parameters(6);
            E3=parameters(7);
            S1=parameters(8);
            S2=parameters(9);
            S3=parameters(10);
            HWHM_0=parameters(11);
            HWHM_T=parameters(12);
            
            fit_func=@(params,cent_energy)obj.lineshape_model(cent_energy,[params(7),params(6),E_fine,T,E1,E2,params(4),params(1),params(2),params(5),params(3),HWHM_T]);
            
            % % %set lower and upper bound
            lb = [0,0,0,0,0,-3000,0];
            ub = [10,10,100,10,0.5,3000,1.5];
            %
            % %set arbitrary intital guess
            p0=[0.2,0.2,2,4,0.9,2050,1];
            %opts=optimoptions(@lsqcurvefit,'StepTolerance',1E-10,'FunctionTolerance',1E-10);
            
            problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',fit_func,'lb',lb,'ub',ub,'xdata',double(center_energy),'ydata',double(center_spectrum));
            ms = MultiStart('PlotFcns',@gsplotbestf);
            [fit_output,error_out] = run(ms,problem,n_multistarts)
            
            fit_params=[fit_output(7),fit_output(6),E_fine,T,E1,E2,fit_output(4),fit_output(1),fit_output(2),fit_output(5),fit_output(3),HWHM_T];
            errormulti=error_out;
            
            %obtain the FWHM for the underlying spectrum.
            FWHM=obj.get_FWHM_model(center_energy,fit_params);
            
            % adding the model parametes to the struct of the specific
            % spectrum.
            [obj.(spectrum_id).(fit_name).fit_params]=fit_params;
            [obj.(spectrum_id).(fit_name).fit_error]=errormulti;
            [obj.(spectrum_id).(fit_name).FWHM]=FWHM;
        end
        
        %% with the ohmic spec dens with exponential cut-off
        function obj=fit_lineshape_prefac_E0_E1_E2_S1_S2_HWHM_exp_spec_dens(obj,spectrum_id,parameters_exp_spec_dens,n_multistarts,fit_name,n,N_satel);
            
            %retrive spectrum from the spectrum_id
            disp(spectrum_id)
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            
            figure()
            plot(energy,spectrum)
            %center and normalize the spectrum
          %  [center_energy,center_spectrum]=center_spec(energy,spectrum);
            center_energy=energy; % not centereing the spectrum because we also want to have the E0;
            center_spectrum=spectrum/max(spectrum);
          
            %% set parameters for fit
            
            prefactor=parameters_exp_spec_dens(1);
            E0=parameters_exp_spec_dens(2);
            E_fine=parameters_exp_spec_dens(3);
            T=parameters_exp_spec_dens(4);
            E1=parameters_exp_spec_dens(5);
            E2=parameters_exp_spec_dens(6);
            cut_off_energy=parameters_exp_spec_dens(7);
            prefac_spec_dens=parameters_exp_spec_dens(8);
            S1=parameters_exp_spec_dens(9);
            S2=parameters_exp_spec_dens(10);
            HWHM_0=parameters_exp_spec_dens(11);
            HWHM_T=parameters_exp_spec_dens(12);
            
            fit_func=@(params,cent_energy) obj.lineshape_model_exp_spec_dens(cent_energy,[params(1),params(2),E_fine,T,params(8),params(9),params(3),params(4),params(5),params(6),params(7),HWHM_T],N_satel);
            
            % % %set lower and upper bound
            lb = [0.7,-3000,0,0,0,0,0,20,33];
            ub = [2,3000,10,5,0.6,0.6,10,30,42];
            %
            % %set arbitrary intital guess
            p0=[1,2010,1,1,0.2,0.2,1,0.2,0.2];
            
            %            opts=optimoptions(@lsqcurvefit,'StepTolerance',1E-10,'FunctionTolerance',1E-10);
            
            problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',fit_func,'lb',lb,'ub',ub,'xdata',double(center_energy),'ydata',double(center_spectrum)');
            ms = MultiStart('PlotFcns',@gsplotbestf);
            [fit_output,error_out] = run(ms,problem,n_multistarts);
            
            fit_params=[fit_output(1),fit_output(2),E_fine,T,fit_output(8),fit_output(9),fit_output(3),fit_output(4),fit_output(5),fit_output(6),fit_output(7),HWHM_T];
            errormulti=error_out;
            
            %obtain the FWHM for the underlying spectrum.
            FWHM=obj.get_FWHM_model(center_energy,fit_params);
            
            % adding the model parametes to the struct of the specific
            % spectrum.
            [obj.(spectrum_id).(fit_name).fit_params]=fit_params;
            [obj.(spectrum_id).(fit_name).fit_error]=errormulti;
            [obj.(spectrum_id).(fit_name).FWHM]=FWHM;
        end
        function obj=fit_lineshape_prefac_E0_HWHM_exp_spec_dens_fixed_HRP(obj,spectrum_id,parameters_exp_spec_dens,n_multistarts,fit_name,n,N_satel);
            
            %retrive spectrum from the spectrum_id
            disp(spectrum_id)
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            
            figure()
            plot(energy,spectrum)
            %center and normalize the spectrum
          %  [center_energy,center_spectrum]=center_spec(energy,spectrum);
            center_energy=energy; % not centereing the spectrum because we also want to have the E0;
            center_spectrum=spectrum/max(spectrum);
          
            %% set parameters for fit
            
            prefactor=parameters_exp_spec_dens(1);
            E0=parameters_exp_spec_dens(2);
            E_fine=parameters_exp_spec_dens(3);
            T=parameters_exp_spec_dens(4);
            E1=parameters_exp_spec_dens(5);
            E2=parameters_exp_spec_dens(6);
            cut_off_energy=parameters_exp_spec_dens(7);
            prefac_spec_dens=parameters_exp_spec_dens(8);
            S1=parameters_exp_spec_dens(9);
            S2=parameters_exp_spec_dens(10);
            HWHM_0=parameters_exp_spec_dens(11);
            HWHM_T=parameters_exp_spec_dens(12);
            
            fit_func=@(params,cent_energy) obj.lineshape_model_exp_spec_dens(cent_energy,[params(1),params(2),E_fine,T,E1,E2,cut_off_energy,prefac_spec_dens,S1,S2,params(3),HWHM_T],N_satel);
            
            % % %set lower and upper bound
            lb = [0,-3000,0];
            ub = [10,3000,100];
            %
            % %set arbitrary intital guess
            p0=[1,1950,1];
            
            %            opts=optimoptions(@lsqcurvefit,'StepTolerance',1E-10,'FunctionTolerance',1E-10);
            
            problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',fit_func,'lb',lb,'ub',ub,'xdata',double(center_energy),'ydata',double(center_spectrum)');
            ms = MultiStart('PlotFcns',@gsplotbestf);
            [fit_output,error_out] = run(ms,problem,n_multistarts);
            
            fit_params=[fit_output(1),fit_output(2),E_fine,T,E1,E2,cut_off_energy,prefac_spec_dens,S1,S2,fit_output(3),HWHM_T];
            errormulti=error_out;
            
            %obtain the FWHM for the underlying spectrum.
            FWHM=obj.get_FWHM_model(center_energy,fit_params);
            
            % adding the model parametes to the struct of the specific
            % spectrum.
            [obj.(spectrum_id).(fit_name).fit_params]=fit_params;
            [obj.(spectrum_id).(fit_name).fit_error]=errormulti;
            [obj.(spectrum_id).(fit_name).FWHM]=FWHM;
        end

        
        %% batch fitting and plotting
        function obj=fit_all_selected_frames(obj,fit_name,parameters,n_multistarts,fit_function_name,n,Nsatel);
            %loop over all structs starting with frame.
            X=0;
            names=fieldnames(obj);
            for i=1:numel(names)
                if strcmp(names(i),'frames')==0
                    if strcmp(names(i),'frameTimes')==0
                        if strcmp(names(i),'framevector')==0
                            id=names(i)
                            if findstr(char(id),'frame')==1
                                X=X+1;
                                if strcmp(fit_function_name,'fit_lineshape_prefac_E0_E1_E2_S1_S2_HWHM_exp_spec_dens')
                                    obj.(fit_function_name)(char(names(i)),parameters,n_multistarts,fit_name,n,Nsatel)
                                elseif strcmp(fit_function_name,'fit_lineshape_prefac_E0_HWHM_exp_spec_dens_fixed_HRP')
                                    obj.(fit_function_name)(char(names(i)),parameters,n_multistarts,fit_name,n,Nsatel)
                                else
                                    obj.(fit_function_name)(char(names(i)),parameters,n_multistarts,fit_name)
                                end
                            end
                        end
                    end
                end
            end
        end
        function obj=plot_all_selected_frames(obj,fit_name,parameters,Nsatel);
            X=0;
            names=fieldnames(obj);
            for i=1:numel(names)
                if strcmp(names(i),'frames')==0
                    if strcmp(names(i),'frameTimes')==0
                        id=names(i)
                        if findstr(char(id),'frame')==1
                            X=X+1;
                        end
                    end
                end
            end
            figure('position', [0, 0, 1400, 900])  % create new figure with specified size
            y=0;
            for i=1:numel(names)
                if strcmp(names(i),'frames')==0
                    if strcmp(names(i),'frameTimes')==0
                        if strcmp(names(i),'framevector')==0
                            id=names(i)
                            if findstr(char(id),'frame')==1
                                y=y+1;
                                axi = subplot(round(X/2)+1,2,y);
                                obj.plot_model_and_lineshape_exp_spec_dens(char(names(i)),fit_name,char(names(i)),Nsatel);
                                ylim([0,1])
                                xlabel('Energy [meV]')
                            end
                        end
                    end
                end
            end
            
        end
        function obj=create_spectra_structs_from_framevector(obj,framevector);
            
            for k=1:numel(framevector);
                spectrum_id=strcat('frame',num2str(framevector(k)));
                obj.create_spectrum_struct(obj.calib_eV,obj.frames(:,framevector(k)),spectrum_id)
            end
            
            if ~isprop(obj,'framevector')% add property to the class
                dummy=addprop(obj,'framevector');
                
            end
            obj.framevector=framevector;
        end
        
        %% plotting.
        function obj=plot_model_and_lineshape(obj,spectrum_id,fit_name,caption)
            
            %retrive spectrum from the spectrum_id
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            model_params=obj.(spectrum_id).(fit_name).fit_params;
            
            %center and normalize the spectrum
            [center_energy,center_spectrum]=center_spec(energy,spectrum);
            center_spectrum=center_spectrum/max(center_spectrum);
            
            lshape=obj.lineshape_model(center_energy,model_params);
            str1=num2str(model_params(1:6),2);
            str2=num2str(model_params(7:12),2);
            str3=caption;
            
            plot(center_energy,center_spectrum);
            hold on
            plot(center_energy,lshape);
            ylim([0,1])
            ylim1=get(gca,'ylim');
            xlim1=get(gca,'xlim');
            disp(ylim1(2)/4)
            %  text(xlim1(1)+xlim1(2)/8,ylim1(2)-ylim1(2)/6,str1,'FontSize',6);
            %  text(xlim1(1)+xlim1(2)/8,ylim1(2)-ylim1(2)/2,str2,'FontSize',6);
            %  text(xlim1(2)-xlim1(2)/2,ylim1(2)-ylim1(2)/2,str3,'FontSize',6);
            
        end
        function obj=plot_model_and_lineshape_exp_spec_dens(obj,spectrum_id,fit_name,caption,N_satel)
            
            %retrive spectrum from the spectrum_id
            energy=obj.(spectrum_id).energy;
            spectrum=obj.(spectrum_id).spectrum;
            model_params=obj.(spectrum_id).(fit_name).fit_params;
            disp(model_params)
            %center and normalize the spectrum
           % [center_energy,center_spectrum]=center_spec(energy,spectrum);
           center_energy=energy; 
           center_spectrum=spectrum/max(spectrum);
            
            lshape=obj.lineshape_model_exp_spec_dens(center_energy,model_params,N_satel);
            str1=num2str(model_params(1:6),2);
            str2=num2str(model_params(7:12),2);
            str3=caption;
            
            plot(center_energy,center_spectrum);
            hold on
            plot(center_energy,lshape);
            ylim([0,1])
            ylim1=get(gca,'ylim');
            xlim1=get(gca,'xlim');
            disp(ylim1(2)/4)
            title(spectrum_id, 'Interpreter','none')
            %  text(xlim1(1)+xlim1(2)/8,ylim1(2)-ylim1(2)/6,str1,'FontSize',6);
            %  text(xlim1(1)+xlim1(2)/8,ylim1(2)-ylim1(2)/2,str2,'FontSize',6);
            %  text(xlim1(2)-xlim1(2)/2,ylim1(2)-ylim1(2)/2,str3,'FontSize',6);
            
        end
        
        %% get average coupling based on fits.
        function obj=get_coupling(obj,framevector,fit_name,model_param_id)
            %framevector is a vector {0:number_of_frames} that contains the
            %lineshapes to be averaged to get the average coupling and the error.
            
            coupling=[];
            names=fieldnames(obj);
            for k=1:numel(framevector);
                str1=strcat('frame',num2str(framevector(k)));
                for i=1:numel(names)
                    if strcmp(names(i),'frames')==0
                        if strcmp(names(i),'frameTimes')==0
                            id=names(i)
                            if findstr(char(id),str1)==1
                                coupling(k,:)=obj.(char(names(i))).(fit_name).fit_params;
                            end
                        end
                    end
                end
            end
            if length(framevector)==1
                params=coupling;
                stderr=zeros(12,1)';
                
            else
                params=mean(coupling);
                stderr=std(coupling);
            end
            if ~isprop(obj,model_param_id)% add property to the class
                dummy=addprop(obj,model_param_id);
            end
            
            obj.(model_param_id)=struct('params',params,'stderr',stderr);
            
            
        end
        
        %% low level stuff
        function [cent_energy,cent_spec]=center_spectrum(obj,energy,frame);
            %find maximum of frame and recenter the spectrum
            [argvalue,argmax]=max(frame);
            cent_energy=energy-energy(argmax);
            cent_spec=frame/max(frame); % just normalize
        end
        function [lineshape]=lineshape_model(obj,E_vector,parameters);
            
            prefactor=parameters(1);
            E0=parameters(2);
            E_fine=parameters(3);
            T=parameters(4);
            E1=parameters(5);
            E2=parameters(6);
            E3=parameters(7);
            S1=parameters(8);
            S2=parameters(9);
            S3=parameters(10);
            HWHM_0=parameters(11);
            HWHM_T=parameters(12);
            
            Kb=8.617E-5;
            lineshape=convolve_modes(E_vector,E0,T,E1,E2,E3,S1,S2,S3,HWHM_0,HWHM_T);
            lineshape=lineshape+convolve_modes(E_vector,E0+E_fine,T,E1,E2,E3,S1,S2,S3,HWHM_0,HWHM_T)*exp(-(E_fine)/1000/(Kb*T));
            lineshape=lineshape'/max(lineshape)*prefactor;
            %plot(E_vector,lineshape/max(lineshape))
            
            
            
            function pop=nq(E,T);
                Kb=8.617E-5;
                pop=(exp(E/1000/(Kb*T))-1).^(-1);
            end
            
            function lorentzian=lorenz(E,E0,HWHM);
                lorentzian=1/pi*(HWHM)./((E-E0).^2+HWHM.^2);
            end
            
            %defining weighting factor for the nth occupation number of a given mode q,
            %and a Huang-Rhys factor S
            
            function Iqn=get_Iqn(Eq,T,gq,n);
                Iqn=((nq(Eq,T)+1)/(nq(Eq,T))).^(n/2)*exp(-(gq*(2*nq(Eq,T)+1)))*besseli(n,2*gq*(nq(Eq,T)*(nq(Eq,T)+1)).^(1/2));%
            end
            
            function spectralDensity=get_spec_dentiy(E);
                spectralDensity=10*np.exp(-0.2*(E-3).^2);
            end
            
            %%convolve all the mdoes
            function intensity=convolve_modes(E_vector,E0,T,E1,E2,E3,S1,S2,S3,HWHM_0,HWHM_T);
                intensity=zeros(length(E_vector),1)';
                HWHM=HWHM_0+HWHM_T*T;
                for n1=-5:1:5;
                    for n2=-5:1:5;
                        for n3=-50:1:50;
                            modeEnergy=n1*E1+n2*E2+n3*E3;
                            intensity=intensity+get_Iqn(E1,T,S1,n1)*get_Iqn(E2,T,S2,n2)*get_Iqn(E3,T,S3,n3)*lorenz(E_vector,E0-modeEnergy,HWHM);
                        end
                    end
                end
            end
        end
        function [FWHM]=get_FWHM_model(obj,center_energy,fit_params)
            
            lshape=obj.lineshape_model(center_energy,fit_params);
            
            [l,o]=max(lshape);
            half_max=l/2;
            
            pp=(lshape-half_max).^2;
            
            for k=1:3;
                [u,t]=min(pp);
                half_indeces(k)=t;
                pp(t)=1E12; %remove the minimum
            end
            
            FWHM=abs(center_energy(min(half_indeces))-center_energy(max(half_indeces)));
        end
        function [lineshape_exp_spec_dens]=lineshape_model_exp_spec_dens(obj,E_vector,parameters_exp_spec_dens,N_satel)
            
            
            prefactor=parameters_exp_spec_dens(1);
            E0=parameters_exp_spec_dens(2);
            E_fine=parameters_exp_spec_dens(3);
            T=parameters_exp_spec_dens(4);
            E1=parameters_exp_spec_dens(5);
            E2=parameters_exp_spec_dens(6);
            cut_off_energy=parameters_exp_spec_dens(7);
            prefac_spec_dens=parameters_exp_spec_dens(8);
            S1=parameters_exp_spec_dens(9);
            S2=parameters_exp_spec_dens(10);
            HWHM_0=parameters_exp_spec_dens(11);
            HWHM_T=parameters_exp_spec_dens(12);
            
            Kb=8.617E-5;
            n=4;
            
            lineshape_exp_spec_dens=convolve_all(E_vector,E0,T,E1,S1,E2,S2,cut_off_energy,prefac_spec_dens,n,HWHM_0,HWHM_T,N_satel);
            % lineshape=lineshape+convolve_all(E_vector,E0+E_fine,T,E1,S1,E2,S2,cut_off_E,prefac_spec_dens,n,HWHM_0,HWHM_T,N_satel)*exp(-(E_fine)/1000/(Kb*T)); %% add the convolved modes of the hiher lying fine-structre state.
            lineshape_exp_spec_dens=lineshape_exp_spec_dens/max(lineshape_exp_spec_dens)*prefactor;
            
            
            function pop=nq(E,T);
                Kb=8.617E-5;
                pop=(exp(E/1000/(Kb*T))-1).^(-1);
            end
            
            function lorentzian=lorenz(E,E0,HWHM);
                lorentzian=1/pi*(HWHM)./((E-E0).^2+HWHM.^2);
            end
            
            %defining weighting factor for the nth occupation number of a given mode q,
            %and a Huang-Rhys factor S
            
            function Iqn=get_Iqn(Eq,T,gq,n);
                Iqn=((nq(Eq,T)+1)/(nq(Eq,T))).^(n/2)*exp(-(gq*(2*nq(Eq,T)+1)))*besseli(n,2*gq*(nq(Eq,T)*(nq(Eq,T)+1)).^(1/2));%
            end
            
            function [energy,exponential_cutoff]=exponential_cutoff_function(drude_energy,prefac,n)
                
                energy=linspace(drude_energy,3*drude_energy,n);
                exponential_cutoff=prefac*energy.*exp(-energy./drude_energy)
                
            end
            
            function lineshape=convolve_all(E_vector,E0,T,E1,S1,E2,S2,cut_off_energy,prefac_spec_dens,n,HWHM_0,HWHM_T,N_satel)
                lineshape=zeros(length(E_vector),1)';
                %N_satel=10;
                [energy,HRP]=exponential_cutoff_function(cut_off_energy,prefac_spec_dens,n)
                %energy=[2,2,2,2];
                %HRP=[0.2,0.2,0.2,0.2];
                HWHM=HWHM_0%+HWHM_T*T;
                for n1=-N_satel:1:N_satel;
                    for n2=-N_satel:1:N_satel;
                        for n3=-N_satel:1:N_satel; %acoustic from spectral density
                            for n4=-N_satel:1:N_satel;%acoustic from spectral density
                                for n5=-3:1:3; %optical 1
                                    for n6=-3:1:3; %optical 2
                                        modeEnergy=n1*energy(1)+n2*energy(2)+n3*energy(3)+n4*energy(4)+n5*E1+n6*E2;
                                        lineshape=lineshape+get_Iqn(energy(1),T,HRP(1),n1)*get_Iqn(energy(2),T,HRP(2),n2)*get_Iqn(energy(3),T,HRP(3),n3)*get_Iqn(energy(4),T,HRP(4),n4)*get_Iqn(E1,T,S1,n5)*get_Iqn(E2,T,S2,n6)*lorenz(E_vector,E0-modeEnergy,HWHM);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end