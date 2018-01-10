%MATLAB Implementation for analysis of photon stream data from the Picoqunt
%Hydraharp.

%% photons.m V3.0 @C HENDRIK UTZAT, KATIE SHULENBERGER, TIMOTHY SINCLAIR 2016
%/16/12/13

%This is an object oriented implementation that allows dynamic
%sorting and analysis of the photon stream collected with the picoquant
%Hydraharp. -NOTE TO MYSELF ADD DESCIPTION-

classdef photons<dynamicprops
    %%%% creating sub-class of dynamicprops which allows
    %to dynamically add propterties to an instance of that class

    %Initialize some properties of photons class.
    properties
        file_path=[];pico_path=[];pathstr=[];fname=[];fext=[];% .ht3 file location specific
        header_size=[]; %size of header in bytes
        header=[];% all header info
        mode=[];%t3 vs t2
        int_time=[]; %integration time
        resolution=[];
        sync_rate=[];
        in_ch0=[];in_ch1=[];in_ch2=[];in_ch3=[]; %average channel counts
        data=[]; %
        n_records=[]; % total number of records in the .ht3 file.
        buffer_size=1E6; %default buffer size for reading in photon records in binary.
        earliest_photon = 0; %for parsing only parts of the photon stream
        latest_photon = 0; %for parsing only parts of the photon stream. 
    end
    %%
    methods %these are build in function that we can perform on the photons object.
        
        %% constructor method
        function obj=photons(file_path,buffer_size,pico_path)
            %%constructor mehtod.
            %It is called automatically when an instance of
            %the photon class is created.It reads the header info from the
            %input .ht3 file. %file_path being the path to the .ht3 file to
            %analyze and need the path to Tom's photon-analysis package.Buffer size
            %is the number of records to load in a batch - typically 1000000.
            
            %set the environmental variable to the path of the photon
            %analysis program of TOM. Will allow to call
            setenv('PATH',pico_path);
            
            %set buffer_size for reading in photons to buffer_size if
            %specified.
            if nargin>1
                obj.buffer_size=buffer_size;
            end
            
            %set path for c_implementation of correlation codes.
            obj.pico_path=pico_path;
            
            %%extract .ht3 file path info.
            [pathstr,name,ext] = fileparts(file_path);
            obj.fname=name;
            obj.file_path=pathstr;
            obj.fext=ext;
            obj.file_path=file_path
            
            %initialize the respective header info for the photon object.
            
            fid=fopen(file_path);
            
            %%code provided from picoquant to read photon arrival times.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % ASCII file header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Ident = char(fread(fid, 16, 'char'));
            fprintf(1,'               Ident: %s\n', Ident);
            
            FormatVersion = deblank(char(fread(fid, 6, 'char')'));
            fprintf(1,'      Format version: %s\n', FormatVersion);
            
            if not(strcmp(FormatVersion,'2.0'))
                fprintf(1,'\n\n      Warning: This program is for version 2.0 only. Aborted.');
                STOP;
            end;
            
            CreatorName = char(fread(fid, 18, 'char'));
            fprintf(1,'        Creator name: %s\n', CreatorName);
            
            CreatorVersion = char(fread(fid, 12, 'char'));
            fprintf(1,'     Creator version: %s\n', CreatorVersion);
            
            FileTime = char(fread(fid, 18, 'char'));
            fprintf(1,'    Time of creation: %s\n', FileTime);
            
            CRLF = char(fread(fid, 2, 'char'));
            
            Comment = char(fread(fid, 256, 'char'));
            fprintf(1,'             Comment: %s\n', Comment);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Binary file header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % The binary file header information is indentical to that in HHD files.
            % Note that some items are not meaningful in the time tagging modes
            % therefore we do not output them.
            
            NumberOfCurves = fread(fid, 1, 'int32');
            
            BitsPerRecord = fread(fid, 1, 'int32');
            fprintf(1,'       Bits / Record: %d\n', BitsPerRecord);
            
            ActiveCurve = fread(fid, 1, 'int32');
            
            MeasurementMode = fread(fid, 1, 'int32');
            fprintf(1,'    Measurement Mode: %d\n', MeasurementMode);
            obj.mode=MeasurementMode
            
            SubMode = fread(fid, 1, 'int32');
            fprintf(1,'            Sub-Mode: %d\n', SubMode);
            
            Binning = fread(fid, 1, 'int32');
            fprintf(1,'             Binning: %d\n', Binning);
            
            Resolution = fread(fid, 1, 'double');
            fprintf(1,'          Resolution: %f ps\n', Resolution);
            obj.resolution=Resolution
            
            Offset = fread(fid, 1, 'int32');
            fprintf(1,'              Offset: %d\n', Offset);
            
            Tacq = fread(fid, 1, 'int32');
            fprintf(1,'    Acquisition Time: %d ms \n', Tacq);
            
            StopAt = fread(fid, 1, 'uint32');
            StopOnOvfl = fread(fid, 1, 'int32');
            Restart = fread(fid, 1, 'int32');
            DispLinLog = fread(fid, 1, 'int32');
            DispTimeAxisFrom = fread(fid, 1, 'int32');
            DispTimeAxisTo = fread(fid, 1, 'int32');
            DispCountAxisFrom = fread(fid, 1, 'int32');
            DispCountAxisTo = fread(fid, 1, 'int32');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i = 1:8
                DispCurveMapTo(i) = fread(fid, 1, 'int32');
                DispCurveShow(i) = fread(fid, 1, 'int32');
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i = 1:3
                ParamStart(i) = fread(fid, 1, 'float');
                ParamStep(i) = fread(fid, 1, 'float');
                ParamEnd(i) = fread(fid, 1, 'float');
            end;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            RepeatMode = fread(fid, 1, 'int32');
            RepeatsPerCurve = fread(fid, 1, 'int32');
            Repaobjime = fread(fid, 1, 'int32');
            RepeatWaiobjime = fread(fid, 1, 'int32');
            ScriptName = char(fread(fid, 20, 'char'));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %          Hardware information header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf(1,'-------------------------------------\n');
            
            HardwareIdent = char(fread(fid, 16, 'char'));
            fprintf(1,' Hardware Identifier: %s\n', HardwareIdent);
            
            HardwarePartNo = char(fread(fid, 8, 'char'));
            fprintf(1,'Hardware Part Number: %s\n', HardwarePartNo);
            
            HardwareSerial = fread(fid, 1, 'int32');
            fprintf(1,'    HW Serial Number: %d\n', HardwareSerial);
            
            nModulesPresent = fread(fid, 1, 'int32');
            fprintf(1,'     Modules present: %d\n', nModulesPresent);
            
            for i=1:10
                ModelCode(i) = fread(fid, 1, 'int32');
                VersionCode(i) = fread(fid, 1, 'int32');
            end;
            for i=1:nModulesPresent
                fprintf(1,'      ModuleInfo[%02d]: %08x %08x\n', i-1, ModelCode(i), VersionCode(i));
            end;
            
            BaseResolution = fread(fid, 1, 'double');
            fprintf(1,'      BaseResolution: %f\n', BaseResolution);
            
            InputsEnabled = fread(fid, 1, 'ubit64');
            fprintf(1,'      Inputs Enabled: %x\n', InputsEnabled); %actually a bitfield
            
            InpChansPresent  = fread(fid, 1, 'int32');
            fprintf(1,' Input Chan. Present: %d\n', InpChansPresent);
            
            RefClockSource  = fread(fid, 1, 'int32');
            fprintf(1,'      RefClockSource: %d\n', RefClockSource);
            
            ExtDevices  = fread(fid, 1, 'int32');
            fprintf(1,'    External Devices: %x\n', ExtDevices); %actually a bitfield
            
            MarkerSeobjings  = fread(fid, 1, 'int32');
            fprintf(1,'     Marker Seobjings: %x\n', MarkerSeobjings); %actually a bitfield
            
            SyncDivider = fread(fid, 1, 'int32');
            fprintf(1,'        Sync divider: %d \n', SyncDivider);
            
            SyncCFDLevel = fread(fid, 1, 'int32');
            fprintf(1,'      Sync CFD Level: %d mV\n', SyncCFDLevel);
            
            SyncCFDZeroCross = fread(fid, 1, 'int32');
            fprintf(1,'  Sync CFD ZeroCross: %d mV\n', SyncCFDZeroCross);
            
            SyncOffset = fread(fid, 1, 'int32');
            fprintf(1,'         Sync Offset: %d\n', SyncOffset);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %          Channels' information header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for i=1:InpChansPresent
                InputModuleIndex(i) = fread(fid, 1, 'int32');
                InputCFDLevel(i) = fread(fid, 1, 'int32');
                InputCFDZeroCross(i) = fread(fid, 1, 'int32');
                InputOffset(i) = fread(fid, 1, 'int32');
                
                fprintf(1,'\n-------------------------------------\n');
                fprintf(1,'Input Channel No. %d\n', i-1);
                fprintf(1,'-------------------------------------\n');
                fprintf(1,'  Input Module Index: %d\n', InputModuleIndex(i));
                fprintf(1,'     Input CFD Level: %d mV\n', InputCFDLevel(i));
                fprintf(1,' Input CFD ZeroCross: %d mV\n', InputCFDZeroCross(i));
                fprintf(1,'        Input Offset: %d\n', InputOffset(i));
            end;
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %                Time tagging mode specific header
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fprintf(1,'\n-------------------------------------\n');
            for i=1:InpChansPresent
                InputRate(i) = fread(fid, 1, 'int32');
                fprintf(1,'     Input Rate [%02d]: %d\n', i-1, InputRate(i));
            end;
            
            fprintf(1,'\-------------------------------------\n');
            
            SyncRate = fread(fid, 1, 'int32');
            fprintf(1,'           Sync Rate: %d Hz\n', SyncRate);
            obj.sync_rate=SyncRate
            
            StopAfter = fread(fid, 1, 'int32');
            fprintf(1,'          Stop After: %d ms \n', StopAfter);
            obj.int_time=StopAfter % given in ms.
            
            StopReason = fread(fid, 1, 'int32');
            fprintf(1,'         Stop Reason: %d\n', StopReason);
            
            ImgHdrSize = fread(fid, 1, 'int32');
            fprintf(1,' Imaging Header Size: %d bytes\n', ImgHdrSize);
            
            nRecords = fread(fid, 1, 'uint64');
            fprintf(1,'   Number of Records: %d\n', nRecords);
            obj.n_records=nRecords;
            
            % Special header for imaging. How many of the following ImgHdr array elements
            % are actually present in the file is indicated by ImgHdrSize above.
            % Storage must be allocated dynamically if ImgHdrSize other than 0 is found.
            
            ImgHdr = fread(fid, ImgHdrSize, 'int32');  % You have to properly interpret ImgHdr if you want to generate an image
            
            % The header section end after ImgHdr. Following in the file are only event records.
            % How many of them actually are in the file is indicated by nRecords in above.
            
            obj.header_size=ftell(fid) %header size is current byte read in .ht3 file
            fclose(fid)
        end
        
        %% normal methods you can call.
        
        %% A) Parsing, chopping and tweaking the primary photon stream data. 
        function obj=picoquant_bin(obj)
            %creates a binary file with all the photon records. The binary file
            %will be accessible via obj.all_photons.photons
            %creates separate binary photons files for each channel used.
            %bugs and to-do: batch write the sorted stream
            %NOTE: It is important to understand the difference between
            %'records' and 'photons'. Each record in the .ht3 file has 32
            %bits encoding the channel, number of sync pulse, the time
            %after the sync and a special character. The special cahracter
            %is important to distinguish whether a record is a true
            %photon. If the special bit is 0, the 32 bit record
            %corresponds to a true photon arrival event. If the record is
            %1, the record is a special overflow record and not a true photon arrival event
            %This means that
            %the 10 bits (1024) were not sufficient to store the number of
            %the sync pulse. By introducing the overflow (OFL) records,
            %the hydraharp software can record the number of syncpulses to
            %infinity, because we can simply account for the number of
            %overflows when reading the .ht3 file to recover the true sync
            %pulse.
            
            %%make distinciton between t2 and t3 data.
            %dealing with t3 (pulsed) data first.
            if obj.mode ==3;
                
                buffer_size=obj.buffer_size;
                header_size=obj.header_size;
                
                
                %for reading the .ht3 binary data file
                syncperiod = 1E9/obj.sync_rate;      % in nanoseconds
                OverflowCorrection = 0;
                T3WRAPAROUND=1024; %if overflow occured, the true n_sync is n_sync+1024
                true_photons=zeros(3,length(buffer_size));%initialize an array to store true_photon records.
                
                outfile = strcat(obj.pathstr,obj.fname,'.photons');
                
                fpout = fopen(outfile,'W');
                fid=fopen(obj.file_path);%opdn the binary .ht3 file
                
                
                fseek(fid,header_size,'bof');%skip over the header to the photon data
                while 1 %while true
                    
                    batch=fread(fid,buffer_size,'ubit32');%reading in a multiple of 32 bit registers
                    lbatch=length(batch);
                    
                    
                    k=0;%true photon counting variable
                    for i=1:lbatch;%looping over all records in batch
                        %read and decode the 32 bit register of the ith record
                        nsync = bitand(batch(i),1023);      % the lowest 10 bits of the ith photon
                        dtime = bitand(bitshift(batch(i),-10),32767);   % the next 15 bits
                        channel = bitand(bitshift(batch(i),-25),63);   % the next 6 bits:%0-4
                        special = bitand(bitshift(batch(i),-31),1);   % the last bit:% MSB - for overflow handling
                        
                        if special == 0   % this means a true 'photon' arrival event.
                            true_nSync = OverflowCorrection + nsync;
                            %  one nsync time unit equals to "syncperiod" which can be calculated from "SyncRate"
                            time =dtime*obj.resolution;
                            k=k+1;%counting the real photons that we see.
                            true_photons(:,k)=[channel;true_nSync;time]; %writing the true photon to an array.
                        else    % this means we have a special record; the 'record' is not a 'photon'
                            if channel == 63  % overflow of nsync occured
                                if(nsync==0) % if nsync is zero it is an old style single oferflow
                                    OverflowCorrection = OverflowCorrection + T3WRAPAROUND;
                                else         % otherwise nsync indicates the number of overflows - THIS IS NEW IN FORMAT V2.0
                                    OverflowCorrection = OverflowCorrection + T3WRAPAROUND*nsync;
                                end;
                            end;
                        end;
                    end;
                    fwrite(fpout,true_photons(:,1:k),'uint64');% writing the true photons to the output file in binary.
                    
                    %break the while loop when we have reached the end of the
                    %.ht3 file.
                    if lbatch <buffer_size;
                        break
                    end
                    
                end
                fclose(fid);
                fclose(fpout);
            end
            
            
            %%here we deal with t2 data.
            if obj.mode==2;
                
                buffer_size=obj.buffer_size;
                header_size=obj.header_size;
                
                %for reading the .ht2 binary data file
                cnt_OFL=0; cnt_MAR=0;  cnt_SYN=0; % just counters
                OverflowCorrection = 0;
                T2WRAPAROUND=33554432; % = 2^25  IMPORTANT! THIS IS NEW IN FORMAT V2.0
                
                true_photons=zeros(2,length(buffer_size));%initialize an array to store true_photon records.
                
                outfile = strcat(obj.pathstr,obj.fname,'.photons');
                
                fpout = fopen(outfile,'W');
                fid=fopen(obj.file_path);%opdn the binary .ht2 file
                
                
                fseek(fid,header_size,'bof');%skip over the header to the photon data
                while 1 %while true
                    
                    batch=fread(fid,buffer_size,'ubit32');%reading in a multiple of 32 bit registers
                    lbatch=length(batch);
                    
                    
                    k=0;%true photon counting variable
                    for i=1:lbatch;%looping over all records in batch
                        
                        %read and decode the 32 bit register of the ith record
                        dtime = bitand(batch(i),33554431);   % the last 25 bits:
                        channel = bitand(bitshift(batch(i),-25),63);   % the next 6 bits:
                        special = bitand(bitshift(batch(i),-31),1);   % the last bit:
                        truetime = OverflowCorrection + dtime;
                        
                        if special == 0   % this means a true 'photon' arrival event.
                            k=k+1;%counting the real photons that we see.
                            true_photons(:,k)=[channel;truetime]; %writing the true photon to a binary array.
                            
                        else    % this means we have a special record; the 'record' is not a 'photon'
                            
                            if channel == 63  % overflow of dtime occured
                                if(dtime==0) % if dtime is zero it is an old style single oferflow
                                    OverflowCorrection = OverflowCorrection + T2WRAPAROUND;
                                    cnt_OFL=cnt_OFL+1;
                                else         % otherwise dtime indicates the number of overflows - THIS IS NEW IN FORMAT V2.0
                                    OverflowCorrection = OverflowCorrection + T2WRAPAROUND*dtime;
                                    cnt_OFL=cnt_OFL+dtime;
                                end;
                            end;
                            
                            if (channel>=1)&(channel<=15);  % these are markers
                                cnt_MAR=cnt_MAR+1;
                                true_photons(:,k)=[99;truetime]; %% if channel = 99, the photon is an external marker.                            end;
                            end;
                        end
                    end
                        fwrite(fpout,true_photons(:,1:k),'uint64');% writing the true photons to the output file in binary.
                        
                        %break the while loop when we have reached the end of the
                        %.ht2 file.
                        if lbatch <buffer_size;
                            break
                        end
                end
                fclose(fid);
                fclose(fpout);
                
            end
        end
        function obj=photons_2_channels(obj,file_in_key,file_out_key,nchannels)
            %% creates four .photons output files containing the photon arrival data of each channel.
            %file_in_key is the filename of the .photons file without
            %extension. file_out_key is the filename (-_ch0 ....3) of the new . photons
            %files.n_channels is the number of channels used (2 for PCFS).
            
            if nargin<4
                nchannels=4;
            end
            
            buffer_size=obj.buffer_size;
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            
            %create output files for each channel
            for i=1:nchannels;
                fid_out(i)=fopen(strcat(obj.pathstr,file_out_key,'_ch_',num2str(i),'.photons'),'w');
            end
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%;'%reading in an array of the binary to process.
                lbatch=length(batch);
                
                for k=1:lbatch;
                    fwrite(fid_out(uint64(batch(k,1))+1),batch(k,:)','ubit64');
                end
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
                
            end
            fclose(fid);
            for i=1:nchannels;
                fclose(fid_out(i));
            end
        end
        function obj=bin_2_int(obj,file_in_key,file_out_key)
            %% converts any picoquant_bin output file into a decimal (human readable), comma separated
            %file for use with Thomas Bischof's correlation code.
            %(https://github.com/tsbischof/photon_correlation)
            
            %file_in_key is the filename of the .photons file without
            %extension. file_out_key is the filename of the created
            %.photons_int
            
            time_min = 0;
            time_max = 0;
            
            buffer_size = obj.buffer_size;
            
            fid = fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            fid_out = fopen(strcat(obj.pathstr,file_out_key,'.photons_int'),'w');
            
            if obj.mode == 3;
            
            foo = 0;
            while 1 %~while true
                foo = foo + 1;
                batch = fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                fprintf(fid_out,'%2u,%9u,%6u\n',batch');
                lbatch = length(batch);
                
                if time_min > min(batch(3,:))
                    time_min = min(batch(3,:));
                end
                if time_max < max(batch(3,:))
                    time_max = max(batch(3,:));
                end
                
                if lbatch < buffer_size; %stop reading in batches when we have reached end of file.
                    break
                end
            end
            fclose(fid)
            fclose(fid_out)
            end
            if obj.mod == 2;
            foo = 0;
            while 1 %~while true
                foo = foo + 1;
                batch = fread(fid,[2,buffer_size],'uint64')';%reading in an array of the binary to process.
                fprintf(fid_out,'%2u,%9u\n',batch');
                lbatch = length(batch);
                
                if time_min > min(batch(2,:))
                    time_min = min(batch(2,:));
                end
                if time_max < max(batch(2,:))
                    time_max = max(batch(2,:));
                end
                
                if lbatch < obj.buffer_size; %stop reading in batches when we have reached end of file.
                    break
                end
            end
            fclose(fid)
            fclose(fid_out)
            obj.earliest_photon = time_min;
            obj.latest_photon = time_max;
            end
            
        end
        function obj=int_file_subset(obj,file_in_key,file_out_key,percent_start,percent_end)
            %% Creates .photons_int files containing only a subset of all photon records
            % it is required that the .bin_2_int() method has already been
            % run so you have a record of the earliest and latest photon.
            % Input is the start and stop times by percent of total time of
            % experiment.            
            % file_in_key is the filename of the .photons file without
            % extension. file_out_key is the filename of the created
            % .photons_int, except the end will be appended with the
            % percent of the total time of experiment that the file's
            % photons represent.
            
            time_start = obj.earliest_photon + (obj.latest_photon - obj.earliest_photon) * percent_start / 100;
            time_end = obj.earliest_photon + (obj.latest_photon - obj.earliest_photon) * percent_end / 100;
            buffer_size = obj.buffer_size;
            fid_in = fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            fid_out = fopen(strcat(obj.pathstr,file_out_key,'_',num2str(percent_start),'_',num2str(percent_end),'.photons_int'),'w');
            if obj.mode == 3
                while 1 %~while true
                    batch = fread(fid_in,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                    lbatch = length(batch);
                    for k = 1:lbatch
                        if (batch(k,3) > time_start) && (batch(k,3) < time_end)
                            fprintf(fid_out,'%2u,%9u,%6u\n',batch(k,:)');
                        end % if
                    end % for
                    if lbatch < buffer_size; %stop reading in batches when we have reached end of file.
                        break
                    end % if
                end % while
            end % if
            if obj.mode == 2
                while 1 %~while true
                    batch = fread(fid_in,[2,buffer_size],'uint64')';%reading in an array of the binary to process.
                    lbatch = length(batch);
                    for k = 1:lbatch
                        if (batch(k,2) > time_start) && (batch(k,2) < time_end)
                            fprintf(fid_out,'%2u,%9u\n',batch(k,:)');
                        end % if
                    end % for
                    if lbatch < obj.buffer_size; %stop reading in batches when we have reached end of file.
                        break
                    end % if
                end % while
                
                fclose(fid_in);
                fclose(fid_out);
            end % if
        end % function
        function obj=write_photons_to_one_channel(obj,file_in_key,file_out_key)
            %writing the photons detected on different channels to one
            %channel. Useful to get the autocorrelation of the sum
            %signal for PCFS analysis.
            %open the binary photon file and change the channel to a
            %new number (0) if the photon was detected at one of the
            %channels given in channels_to_combine.
            if obj.mode==3;
            buffer_size=obj.buffer_size;
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            
            %create output file for the combined signal
            fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w');
            
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%;'%reading in an array of the binary to process.
                lbatch=length(batch);
                batch(:,1)=zeros(lbatch,1);
                fwrite(fid_out,batch','ubit64');
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
            end
            fclose(fid);
            fclose(fid_out);
            end
            
            %dealing with t2 data.
            if obj.mode == 2;
            
            buffer_size=obj.buffer_size;
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            
            %create output file for the combined signal
            fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w');
            
            while 1 %~while true
                batch=fread(fid,[2,buffer_size],'uint64')';%;'%reading in an array of the binary to process.
                lbatch=length(batch);
                batch(:,1)=zeros(lbatch,1);
                fwrite(fid_out,batch','ubit64');
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
            end
            fclose(fid);
            fclose(fid_out);            
            end
        end
        
        %% B) analysis of the photon stream, lifetimes, correlations etc.
        function obj=get_intensity_w_SPEmarkers(obj,file_in_key,file_out_key,bin_width,upper_limit,lower_limit)
            %stores the intensity trace as a property of the photons class.
            %(obj.file_out_key) file_in_key is the filename of the .photons
            %file that is used to compile the intensity trace from. 
            
            %If upper, lower_limit and file_out_key are given, the photons
            %emitted in bins  with intensities in these bounds are written
            %to a file_out_key.photons.
            
            %bin_width is ps for t2 data and number of pulses for t3 data.
            
            %currently only T2 data.
            %%dealing with t2 data
            
            
            
            if obj.mode==3;
                if nargin<5
                    photon_sorting=false;
                else
                    photon_sorting=true;
                end
                
                if photon_sorting==true
                    
                    buffer_size=obj.buffer_size;
                    
                    %counters.
                    bin_count=0;
                    run_count=[0 0 0 0];%counting variable for all four channels.
                    counts=[];
                    
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w');
                    
                    while 1 %~while true
                        batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary photons to process.
                        lbatch=length(batch);
                        indices=[]; %vector storing the indices of photons in a particular bin.
                        for i=1:lbatch;
                            while 1 % ~ while true
                                bin=[bin_count,bin_count+1]*bin_width; %define a new bin in units of pulses
                                if batch(i,2)<= bin(2) & batch(i,2) > bin(1);
                                    run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                    indices=[indices;i];%store photons indeces that are in the bin.
                                    break
                                else %if the photon is not in the bin in units of pulses
                                    bin_count=bin_count+1; %move to the next bin
                                    counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                    if sum(run_count)<=upper_limit & sum(run_count)>=lower_limit & ~isempty(indices)
                                        buffer=batch(indices(1):indices(end),:);%buffer the photons from the last bin
                                        fwrite(fid_out,buffer','ubit64'); %buffer needs to be transposed to match dim of records variable in picoquant_bin.
                                    end
                                    run_count=[0 0 0 0];%reset the counters for the next bin.
                                    indices=[];
                                end
                              
                                
                            end
                        end
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break
                        end
                    end
                    fclose(fid)
                    fclose(fid_out)
                    
                    %store the counts as property of the photons object
                    if ~isprop(obj,file_out_key)
                        dummy=addprop(obj,file_out_key);
                    end
                    obj.(file_out_key)=struct('trace',counts,'bin_width',bin_width,'upper_limit',upper_limit,'lower_limit',lower_limit);
                    
                else
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    buffer_size=obj.buffer_size;
                    
                    bin_count=0;
                    run_count=[0 0 0 0];
                    counts=[];
                    
                    while 1 %~while true
                        batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                        lbatch=length(batch);
                        for i=1:lbatch;
                            while 1 % ~ while true
                                window=[bin_count,bin_count+1]*bin_width; %define a new window in units of pulses
                                if batch(i,2)<= window(2) & batch(i,2) > window(1);
                                    run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                    break
                                else %if the photon is not in the window in units of pulses
                                    bin_count=bin_count+1; %move to the next bin
                                    counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                    run_count=[0 0 0 0];%reset the running count to zero; reiterate until the ith photon is fit into a window
                                end
                            end
                        end
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break
                        end
                    end
                    fclose(fid)
                    if ~isprop(obj,file_out_key)
                        dummy=addprop(obj,file_out_key);
                    end
                    obj.(file_out_key)=struct('trace',counts,'bin_width',bin_width,'upper_limit','n.a','lower_limit','n.a');
                end
                time=linspace(0,length(obj.(file_out_key).trace),length(obj.(file_out_key).trace))*bin_width/obj.sync_rate;
            end
            
            
            if obj.mode==2;
                
                if nargin<5
                    photon_sorting=false;
                else
                    photon_sorting=true;
                end
                
                if photon_sorting==true
                    
                    buffer_size=obj.buffer_size;
                    
                    %counters.
                    bin_count=0;
                    run_count=[0 0 0 0];%counting variable for all four channels.
                    counts=[];
                    
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w');
                    
                    while 1 %~while true
                        batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary photons to process.
                        lbatch=length(batch);
                        indices=[]; %vector storing the indices of photons in a particular bin.
                        for i=1:lbatch;
                            while 1 % ~ while true
                                bin=[bin_count,bin_count+1]*bin_width; %define a new bin in units of pulses
                                if batch(i,2)<= bin(2) & batch(i,2) > bin(1);
                                    run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                    indices=[indices;i];%store photons indeces that are in the bin.
                                    break
                                else %if the photon is not in the bin in units of pulses
                                    bin_count=bin_count+1; %move to the next bin
                                    counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                    if sum(run_count)<=upper_limit & sum(run_count)>=lower_limit & ~isempty(indices)
                                        buffer=batch(indices(1):indices(end),:);%buffer the photons from the last bin
                                        fwrite(fid_out,buffer','ubit64'); %buffer needs to be transposed to match dim of records variable in picoquant_bin.
                                    end
                                    run_count=[0 0 0 0];%reset the counters for the next bin.
                                    indices=[];
                                end
                                
                                
                            end
                        end
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break
                        end
                    end
                    fclose(fid)
                    fclose(fid_out)
                    
                    %store the counts as property of the photons object
                    if ~isprop(obj,file_out_key)
                        dummy=addprop(obj,file_out_key);
                    end
                    obj.(file_out_key)=struct('trace',counts,'bin_width',bin_width,'upper_limit',upper_limit,'lower_limit',lower_limit);
                
                
                else %% we don't sort photons
                    
                    fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
                    buffer_size=obj.buffer_size;
                    
                    bin_count=0;
                    run_count=[0 0 0 0];
                    run_count_markers=0;
                    counts=[];
                    marker_times=[];
                    intensity_between_markers=[];%get integrated intensity between markers. 
                    marker_count=0;
                    
                    while 1 %~while true
                        batch=fread(fid,[2,buffer_size],'uint64')';%reading in an array of the binary to process.
                        lbatch=length(batch);
                        for i=1:lbatch;
                            if batch(i,1)~=99 %% check if the record is a marker. 
                            run_count_markers=run_count_markers+1;
                            while 1 % ~ while true
                                window=[bin_count,bin_count+1]*bin_width; %define a new window in units of ps
                                if batch(i,2)<= window(2) & batch(i,2) > window(1);
                                    run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                    break
                                else %if the photon is not in the window in units of ps
                                    bin_count=bin_count+1; %move to the next bin
                                    counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                    run_count=[0 0 0 0];%reset the running count to zero; reiterate until the ith photon is fit into a window
                                end
                            end
                            else
                                marker_count=marker_count+1;
                                intensity_between_markers(marker_count)=run_count_markers; %add counts between markers
                                run_count_markers=0; %reset marker. 
                                marker_times=[marker_times;batch(i,2)];
                            end
                        end
                        if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                            break
                        end
                    end
                    fclose(fid)
                    if ~isprop(obj,file_out_key)
                        dummy=addprop(obj,file_out_key);
                    end
                    marker_times=marker_times/1E12;
                    obj.(file_out_key)=struct('trace',counts,'bin_width',bin_width,'marker_times',marker_times,'marker_trace',intensity_between_markers,'upper_limit','n.a','lower_limit','n.a');
                end
                
            end
            
            % get the time vector for the binned intensity trace
            time=linspace(0,length(obj.(file_out_key).trace),length(obj.(file_out_key).trace))*bin_width/1E12
            if ~isprop(obj,'time')
                dummy=addprop(obj,'time');
            end
            obj.time=time;
            
            %%plot the intensity trace with units of time
%             figure()
%             plot(obj.time,obj.(file_out_key).trace)
%             xlabel('Time [ps]')
%             ylabel('counts per bin')
%             title(strcat(obj.fname,' - Intensity Trace'))
            
                end
        function rmv_AP(obj,file_in_key,file_out_key,num_det)

            
            buffer_size=obj.buffer_size; %number of photons processed per batch
            
                      
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the photon stream. this should not be sorted by channel
            fid_out=fopen(strcat(obj.pathstr,file_out_key,'.photons'),'w'); %open the file which we write the new photons to
            u = 0;
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                u = u+1; %counter for batches
                %batch_overflow(4:6,:) = batch(1:3,:);
                
                lbatch=length(batch);
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
                skip_next = 0; %originate counting values
                skip = [1:num_det];
                

                for k=1:lbatch-3; % this loop is the crux of the photon sorting. exmaine each pulse, decide how many photons arrive in a pulse, then write only the first photon per channel per pulse
                      
                    if skip_next > 0 %so we don't double count photons. This will allow us to go pulse by pulse
                        skip_next = skip_next-1;
                        continue
                    end
 
                                        
                    while batch(k,2) == batch(k+same_pulse,2) % decides if its still the same pulse
                        pulse_photons(same_pulse,:) = batch(k,:); %checks pulse
                        skip_next = skip_next+1; % determines the photon number event
                        same_pulse = same_pulse+1;
                    end
                    
                    same_pulse = 1; % reset same_pulse for next pulse
                    
                    [a,b] = size(pulse_photons);
                    
                    if a <= 1
                        buffer(count,:) = pulse_photons(find_first,:); %write this first photon to a buffer
                        count = count+1; % adds to the next entry for buffer
                    else
                        for ii = 1:num_det %loop over all possible detector numbers
                            find_first = find(pulse_photons(:,1)==ii,1); % find the first photon per detector within this pulse
                            %write photon from this index of pulse_photons to
                            %the final string
                            buffer(count,:) = pulse_photons(find_first,:); %write this first photon to a buffer
                            count = count+1; % adds to the next entry for buffer

                            clear find_first %reset the index for the first photon per detector
                        end
                    end
                  
                    fwrite(fid_out,buffer','ubit64'); %buffer needs to be transposed to match dim of records variable in picoquant_bin.
                    
                end

            
            end
            fclose(fid)
            fclose(fid_out)

        end
        function obj=FLID(obj,file_in_key,file_out_key,bin_divisor,bin_width_photons,fit_it)
            %BUGLIST: Number of bins set in hh software should be allowed
            %that are not 2^15.
            if obj.mode==2;
                fprint('This function is not available for t2 data.')
            end
            %returns intensity trace, and lifetimes for a binned or the
            %entire photon stream. Output is returned to struct array.
            %optionally, the lifetime data can be fit with an
            %exponential - lifetime parameter is returned in
            %~.fit_params.
            %with obj.file_out_key.counts ~.lifetimes ~.fit_params
            
            %read in large batch of photons.-
            %bin them into bundles with max span of bin_width_photons.
            %histogram them with bin-width defined,
            %make sure last batch bin_width is the same as for all others.
            %add all histograms of all batches which will be lifetime.
            %create propoer time x-axis by taking resolution into account.
            
            buffer_size=obj.buffer_size;
            
            if nargin<5 %
                bin_width_photons=1;
            end
            
            if nargin<6 % by default, we do not fit the lifetimes for each bin.
                fit_it=false;
            end
            
            %counters
            bin_count=0;
            run_count=[0 0 0 0];
            
            %initializations
            histo=[]; % histogramming array for lifetime compilation.
%             true_nbins=2^15/bin_divisor;%the true number of bins is the systems bin depth divided by an integer defined in bin_divisor to avoid artifacts.
%             bin_vector=linspace(0.5,true_nbins,(true_nbins-0.5))*obj.resolution;% center positions for bins.

            rep_rate = 1/obj.sync_rate*10^15;
            bins = rep_rate/(obj.resolution);
            bin_vector = linspace(obj.resolution/2,rep_rate-obj.resolution/2,bins);
            time=bin_vector;%
            histo = zeros(size(bin_vector));
%             max_time = 1/obj.sync_rate*(10^15);
%             time=linspace(0,max_time,2^15);
%             bin_vector = time;
            fit_param=[]; %array to store fit_params.
            
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            
            if bin_width_photons ~=1 %we want to get the traces for a number of bins.
                while 1
                    batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                    
                    lbatch=length(batch);
                    
                    
                    indices=[];%count vector to convert photons into a buffer.
                    for i=1:lbatch; %looping over all photons in batch.
                        while 1 % ~ while true
                            bin=[bin_count,bin_count+1]*bin_width_photons; %define a new window in units of pulses
                            
                            if batch(i,2)<= bin(2) & batch(i,2) > bin(1);
                                run_count(batch(i,1)+1)=run_count(batch(i,1)+1)+1;% add the photon
                                indices=[indices;i];%store photons indeces that are in the bin.
                                break
                            else %if the photon is not in the window in units of pulses
                                bin_count=bin_count+1; %move to the next bin
                                counts(bin_count,:)=run_count; %add the counts of the previous bin to the counts variable
                                if ~isempty(indices)
                                    buffer=batch(indices(1):indices(end),:);
                                    histo(:,bin_count)=hist(buffer(:,3),bin_vector);
                                    if fit_it==true
                                        x=linspace(0,true_nbins,true_nbins)';
                                        [dummy,index]=max(histo(:,bin_count));
                                        f=fit(x(index:end),histo(index:end,bin_count),'Exp1');
                                        fit_param=[fit_param;f.b];
                                    end
                                    %buffer the photons from the last bin
                                else
                                    histo(:,bin_count)=zeros(1,true_nbins);
                                end
                                run_count=[0 0 0 0];%reset counter.
                                indices=[];%reset index vector of photon buffer.
                            end
                        end
                    end
                    if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                        break
                    end
                end
                
                
            else %we just compile a histogram of all photons without temporal binning.
                u=1;%counter for batch
                while 1
                    batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                    u=u+1;%bin counter need to change this for bin_count as previously.
                    lbatch=length(batch);
                    
                    histo=histo+hist(batch(:,3),bin_vector);
                    if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                        break
                    end
                end
                %histo=sum(run_histo,2);%umming up the histograms
                
                %we don't get an intensity trace or fit parameters now.
                counts='n.a';
                fit_param='n.a';
            end
            fclose(fid)
            
            time = time(1:bins/1000);
            histo = histo(1:bins/1000);
            
            %store returns as properties of the photons object.
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('time',time,'lifetime',histo,'intensity',counts,'fit_param',fit_param);
        end
        function obj=FLIM(obj,file_in_key,file_out_key)
            intensity = sum(obj.CdSeCS_NPL_1e4_dot_005_240nW_000_int.trace,2);
            l_int = size(intensity);
            rep_rate = 1/obj.sync_rate*10^12;
            bins = rep_rate/(obj.resolution);
            bin_vector = linspace(obj.resolution/2,rep_rate-obj.resolution/2,bins);
            %time=bin_vector;%
            max_int = max(intensity);
            n_int_bins = 73;
            int_bins = linspace(0,max_int,n_int_bins);
            FLIM_data = zeros(uint16(bins),uint16(n_int_bins));
            %size(FLIM_data)
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r');
            bin_norm = zeros(1,uint16(n_int_bins));
            int = linspace(1,n_int_bins,n_int_bins)*max_int/n_int_bins;

            for ii = 1:l_int %~while true
                clear run_histo_temp
                batch_size = intensity(ii);
                batch=fread(fid,[3,batch_size],'uint64')';%reading in an array of the binary to process.
                lbatch=length(batch);

                run_histo_temp = hist(batch(:,3),bin_vector)';
                %size(run_histo_temp)
                if lbatch <batch_size; %stop reading in batches when we have reached end of .ht3 file.
                    break
                end
                %find way to select int bin for specific intensity
                index = find(int_bins>=intensity(ii),1);
                bin_norm(index) = bin_norm(index)+1;
                %write temporary hist to that trace 
                %size(FLIM_data(:,index))
                %size(run_histo_temp)
                %ii
                %intensity(ii)
                FLIM_data(:,index) = FLIM_data(:,index)+run_histo_temp;
            end
            
            for jj = 1:n_int_bins
                if bin_norm(jj) == 0
                    continue
                end
                FLIM_norm(:,jj) = FLIM_data(:,jj)./bin_norm(jj);
                [max_val, max_ind] = max(FLIM_norm(:,jj));
                tau_half_val(jj) = max_val/50;
                index_th(jj) = find(FLIM_norm(max_ind:end,jj)<=max_val./50,1);
                tau_half(jj) = bin_vector(index_th(jj));
            end
            
            figure;
            plot(int,tau_half/1000)
            figure;
            plot(index_th)
            figure;
            plot(tau_half_val)
            figure;
            plot(bin_vector)
            
            figure;
            plot(bin_norm)
            Heat_Object = HeatMap(FLIM_norm','Colormap',hsv,'ColumnLabels',bin_vector,'RowLabels',int_bins);
            addYLabel(Heat_Object,'Intensity')
            addXLabel(Heat_Object,'\tau (ns)')
            %colormap winter
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('FLIM',FLIM_norm,'Bin_Counts',bin_norm,'Intensity_bins',int);
        end
        function obj=get_g2(obj,file_in_key,file_out_key,n_pulses,time_bins,rrate)
            %rrate = obj.sync_rate;
            %time_bins = 1/obj.sync_rate;
            time_bins = time_bins*n_pulses;

            buffer_size=obj.buffer_size; %number of photons processed per batch
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the phooton stream. this should not be sorted by channel
            u = 0;
            %g2 = zeros(2*n_pulses+1,channels,channels,time_bins*2+1); 
           
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                u = u+1; %counter for batches
                %batch_overflow(4:6,:) = batch(1:3,:);
                
                lbatch=length(batch);
                clear g2_temp
                
                ii = 0;

                while 1 
                    ii = ii +1;
                    batch_shift_down = batch;
                    batch_shift_down(1:ii,:) = [];
                    %size(batch_shift_down)
                    batch_shift_up = batch;
                    batch_shift_up(end-ii+1:end,:) = [];
                    %size(batch_shift_up)
                    batch_diff = batch_shift_down-batch_shift_up;
                    pulse_sorted = batch_diff;
                    batch_diff_extended = [batch_diff,batch_shift_down(:,3),batch_shift_up(:,3)];
                    Logic1 = batch_diff(:,2) > n_pulses; %If the pulse separation is greater than the max
                    Logic2 = batch_diff(:,2) < -n_pulses; %if the pulse separation is less than the neg of the max

                    Logic = Logic1 | Logic2; %Add up all the conditions (could easily add more). This just cuts off pulse separation
                    pulse_sorted(Logic,:) = []; %Remove those rows
                    pulse_sorted_AP = pulse_sorted;
                    Logic3 = pulse_sorted(:,1) == 0; %If they arrive at the same detector
%                     Logic4 = pulse_sorted(:,2) == 0; %Same pulse
                    Logic_AP = Logic3; %& Logic4;
                    pulse_sorted_AP(Logic_AP,:) = []; %Removes After Pulsing
                    
                    %pulse_sorted is now a matrix containing difference in
                    %detector number, difference in time, and pulse
                    %difference. Can now turn pulse difference and time
                    %difference into total time separation and histogram
                    %and plot.
                    pulse_sorted_time = pulse_sorted_AP(:,3)+1/rrate*10^12*(pulse_sorted_AP(:,2));
                    g2_temp(:,ii) = hist(pulse_sorted_time,time_bins*(2*n_pulses+1));
                    if isempty(pulse_sorted_AP)
                        break;
                    end
                    
                                        

                end
                %size(g2_temp)
                g2_batch(u,:) = sum(g2_temp,2);
                
                
                
                
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end                    
            end
            g2 = sum(g2_batch,1);
            size(g2)
            time = (0:(n_pulses+1)*10^9/(rrate*time_bins*(2*n_pulses+1)):(n_pulses+1)*(10^9/rrate)-10^9/(rrate*time_bins*(2*n_pulses+1)));
            size(time)
            figure;
            plot(time(2:end),g2(2:end),'b')
            hold on
            plot(-time,fliplr(g2'),'b')
            center = 2*trapz(g2(1:end/4));
            side = trapz(g2(end/4:end));
            ratio = center/side
            noise = g2(7*end/8:8*end/8);
            average_noise = mean(noise);
            g2_noise_sub = g2-average_noise;
            figure;
            plot(time,g2_noise_sub)
            center_noise = 2*trapz(g2_noise_sub(1:end/4));
            side_noise = trapz(g2_noise_sub(end/4:end));
            ratio_noise = center_noise/side_noise
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('time',time,'g2',g2)
            

        end
        function get_PNRL(obj,file_in_key,file_out_key,bin_divisor)
            % sort photon stream into PNRL streams for 1-4 photons per
            % pulse. This can then be plotted as a lifetime where each of
            % the "channels" corresponds to a photon number and order.
            
            buffer_size=obj.buffer_size; %number of photons processed per batch
            
            %LIFETIME

            %initializations
            true_nbins=2^15/bin_divisor;%the true number of bins is the systems bin depth divided by an integer defined in bin_divisor to avoid artifacts.
            bin_vector=linspace(0.5,true_nbins,(true_nbins-0.5))*obj.resolution*bin_divisor;% center positions for bins.
            time=bin_vector*10^(-3);
            
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the phooton stream. this should not be sorted by channel
            u = 0;
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                u = u+1; %counter for batches
                %batch_overflow(4:6,:) = batch(1:3,:);
                
                lbatch=length(batch);
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end
                skip_3 = 0; %originate counting values
                skip_2 = 0;
                skip_1 = 0;
                count4 = 0;
                count3 = 0;
                count2 = 0;
                PNRL_4_4_temp = [0 0 0];%create temporary arrays
                PNRL_4_3_temp = [0 0 0];
                PNRL_4_2_temp = [0 0 0];
                PNRL_4_1_temp = [0 0 0];
                PNRL_3_3_temp = [0 0 0];
                PNRL_3_2_temp = [0 0 0];
                PNRL_3_1_temp = [0 0 0];
                PNRL_2_2_temp = [0 0 0];
                PNRL_2_1_temp = [0 0 0];
                

                for k=1:lbatch-3; % this loop is the crux of the photon sorting. exmaine each pulse, decide how many photons arrive in a pulse, then write to channels
                    % we need n(n+1)/2 channels to describe the number of n
                    % photon events possible. The default for my setup is
                    % 4. which is what I am writing this for.
                    %read in a photon.
                    %determine pulse number
                    %read next photon
                    %determine pulse number
                    % if same pulse store both and repeat
                    %continue until photon is on later pulse
                    %determine how many total photons were detected. Write
                    %to appropriate PNRL channel
                    %counters
                    
                    %Continue loop if we have identified a MX event so that
                    %we have no overlap
                   
                    if skip_3 > 0 %for 4 photon events so we don't double count
                        skip_3 = skip_3-1;
                        continue
                    end
                    if skip_2 > 0 %for 3 photon events so we don't double count
                        skip_2 = skip_2-1;
                        continue
                    end
                    if skip_1 > 0 % for 2 photon events so we don't double count
                        skip_1 = skip_1-1;
                        continue
                    end
                        
                    
                    %PNR SORTING

                    n_sync = batch(k,2); %determine the pulse number of the photon in question

                
                    if n_sync == batch(k+3,2) %if the fourth next photon was from the same pulse (4 photon event)
                        %write this photon and the previous 3 to PNRL
                        %channels 1-4
                        batch_hist = hist(batch(k:k+3,1),4);
                        if max(batch_hist) > 1%remove after pulsing
                            clear batch_hist
                            continue
                        end
                        count4 = count4+1;
                        PNRL_4_4_temp(count4,:) = batch(k,:);
                        PNRL_4_3_temp(count4,:) = batch(k+1,:);
                        PNRL_4_2_temp(count4,:) = batch(k+2,:);
                        PNRL_4_1_temp(count4,:) = batch(k+3,:);
                        %disp('4 Photon')
                        skip_3 = 3; %to avoid multicounting
                    elseif n_sync == batch(k+2,2) %three photon event
                        %write this photon and the previous 2 to PNRL
                        %channels 5-7
                        batch_hist = hist(batch(k:k+2,1),4);
                        if max(batch_hist) > 1%remove after pulsing
                            clear batch_hist
                            continue
                        end
                        count3 = count3+1;
                        PNRL_3_3_temp(count3,:) = batch(k,:);
                        PNRL_3_2_temp(count3,:) = batch(k+1,:);
                        PNRL_3_1_temp(count3,:) = batch(k+2,:);
                        %disp('3 Photon')
                        skip_2 = 2;
                    elseif n_sync==batch(k+1,2) %two photon event
                        %write this photon and the previous to PNRL
                        %channels 8-9
                        if batch(k,1) == batch(k+1,1)
                            continue
                        end
                        count2=count2+1;
                        PNRL_2_2_temp(count2,:) = batch(k,:);
                        PNRL_2_1_temp(count2,:) = batch(k+1,:);
%                         disp('2 Photon')
%                         disp(count2)
                        skip_1 = 1;
                    end
                    
                    %store returns as properties of the photons object.
%                     if ~isprop(obj,file_out_key)
%                         dummy=addprop(obj,file_out_key)
%                     end
                    
                end
                run_histo(:,1,u)=hist(PNRL_4_4_temp(:,3),bin_vector);
                run_histo(:,2,u)=hist(PNRL_4_3_temp(:,3),bin_vector);
                run_histo(:,3,u)=hist(PNRL_4_2_temp(:,3),bin_vector);
                run_histo(:,4,u)=hist(PNRL_4_1_temp(:,3),bin_vector);
                run_histo(:,5,u)=hist(PNRL_3_3_temp(:,3),bin_vector);
                run_histo(:,6,u)=hist(PNRL_3_2_temp(:,3),bin_vector);
                run_histo(:,7,u)=hist(PNRL_3_1_temp(:,3),bin_vector);
                run_histo(:,8,u)=hist(PNRL_2_2_temp(:,3),bin_vector);
                run_histo(:,9,u)=hist(PNRL_2_1_temp(:,3),bin_vector);
                clear count2
                clear count3
                clear count4
                clear PNRL_4_4_temp
                clear PNRL_4_3_temp
                clear PNRL_4_2_temp
                clear PNRL_4_1_temp
                clear PNRL_3_3_temp
                clear PNRL_3_2_temp
                clear PNRL_3_1_temp
                clear PNRL_2_2_temp
                clear PNRL_2_1_temp
                %clear batch_overflow
                %batch_overflow(1:3,:) = batch(lbatch-2:lbatch,:);               
            end
            [a,b] = size(run_histo);
            PNRL = zeros(a,b);
            PNRL(:,1)=sum(run_histo(:,1,:),3);
            PNRL(:,2)=sum(run_histo(:,2,:),3);
            PNRL(:,3)=sum(run_histo(:,3,:),3);
            PNRL(:,4)=sum(run_histo(:,4,:),3);
            PNRL(:,5)=sum(run_histo(:,5,:),3);
            PNRL(:,6)=sum(run_histo(:,6,:),3);
            PNRL(:,7)=sum(run_histo(:,7,:),3);
            PNRL(:,8)=sum(run_histo(:,8,:),3);
            PNRL(:,9)=sum(run_histo(:,9,:),3);
            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key);
            end
            obj.(file_out_key)=struct('time',time,'PNRL',PNRL);

        end
        function obj=get_PNRL2(obj,file_in_key,file_out_key)
            test_counter = 0;

            buffer_size=obj.buffer_size; %number of photons processed per batch
            fid=fopen(strcat(obj.pathstr,file_in_key,'.photons'),'r'); % open the phooton stream. this should not be sorted by channel
            u = 0;
            %g2 = zeros(2*n_pulses+1,channels,channels,time_bins*2+1); 
            
%             true_nbins=2^15/bin_divisor;%the true number of bins is the systems bin depth divided by an integer defined in bin_divisor to avoid artifacts.
%             bin_vector=linspace(0.5,true_nbins,(true_nbins-0.5))*obj.resolution*bin_divisor;% center positions for bins.
%             max_time = 1/obj.sync_rate*(10^15);
%             time=linspace(0,max_time,2^15);
            rep_rate = 1/obj.sync_rate*10^15;
            bins = rep_rate/(obj.resolution);
            bin_vector = linspace(obj.resolution/2,rep_rate-obj.resolution/2,bins);
            %true_nbins=2^15;%/bin_divisor;%the true number of bins is the systems bin depth divided by an integer defined in bin_divisor to avoid artifacts.
            %bin_vector=linspace(0.5,true_nbins,(true_nbins-0.5))*obj.resolution;% center positions for bins.
            time=bin_vector;%
            size(time)
            figure;
            plot(time)
            PNRL21 = zeros(size(time));
            PNRL22 = zeros(size(time));
            PNWT22 = zeros(size(time));
           
            while 1 %~while true
                batch=fread(fid,[3,buffer_size],'uint64')';%reading in an array of the binary to process.
                u = u+1; %counter for batches
                %batch_overflow(4:6,:) = batch(1:3,:);
                
                lbatch=length(batch);
                clear g2_temp
                
                ii = 0;
                photons_temp = [];
                photons_ext_temp = [];

                while 1 
                    clear pulse_sorted_AP
                    clear batch_diff
                    ii = ii +1;
                    batch_shift_down = batch;
                    batch_shift_down(1:ii,:) = [];
                    %size(batch_shift_down)
                    batch_shift_up = batch;
                    batch_shift_up(end-ii+1:end,:) = [];
                    %size(batch_shift_up)
                    batch_diff = batch_shift_down-batch_shift_up;
                    pulse_sorted = batch_diff;
                    batch_diff_extended = [batch_diff,batch_shift_down(:,2:3),batch_shift_up(:,3)];
                    batch_diff_extended_t1 = batch_diff_extended;
%                     size(batch_diff_extended)
                    Logic1 = batch_diff(:,2) ~= 0; %If the photons don't arrive on the same pulse

                    Logic = Logic1; %Add up all the conditions (could easily add more). This just cuts off pulse separation
                    pulse_sorted(Logic,:) = []; %Remove those rows so we only have multiple photons from the same pulse
                    batch_diff_extended_t1(Logic,:) = [];
                    pulse_sorted_AP = pulse_sorted;
                    Logic3 = pulse_sorted(:,1) == 0; %If they arrive at the same detector
%                     Logic4 = pulse_sorted(:,2) == 0; %Same pulse
                    Logic_AP = Logic3; %& Logic4;
                    pulse_sorted_AP(Logic_AP,:) = []; %Removes After Pulsing
                    batch_diff_extended_t1(Logic_AP,:) = [];
                    
                    %pulse_sorted is now a matrix containing difference in
                    %detector number, difference in time, and pulse
                    %difference. Can now turn pulse difference and time
                    %difference into total time separation and histogram
                    %and plot.
                    photons_ext_temp = [photons_ext_temp;batch_diff_extended_t1];
                    size(photons_ext_temp);
                    photons_temp = [photons_temp;pulse_sorted_AP];
                    if isempty(pulse_sorted_AP)
                        break;
                    end
                                                   
                end
                unique_pulse = unique(photons_temp);%determine which pulses had multiple photons
                unique_pulse_ext = unique(photons_ext_temp(:,4));
                test_counter = test_counter+ size(unique_pulse_ext);
                [count_pulses dummy] = size(unique_pulse_ext); %how many unique pulses
                PNRL21_batch = [];
                PNRL22_batch = [];
                PNWT22_batch = [];
                for jj = 1:count_pulses
                    clear num_ph_temp
                    PNRL21_temp = [];
                    PNRL22_temp = [];
                    PNWT22_temp = [];
                    num_ph_temp = nnz(photons_ext_temp(:,4) == unique_pulse_ext(jj));
                    if num_ph_temp == 1
                        index_temp = find(photons_ext_temp(:,4) == unique_pulse_ext(jj),1);
                        PNRL21_temp = photons_ext_temp(index_temp,6);
                        PNRL22_temp = photons_ext_temp(index_temp,5);
                        PNWT22_temp = photons_ext_temp(index_temp,3);
                    end
                    PNRL21_batch = [PNRL21_batch;PNRL21_temp];
                    PNRL22_batch = [PNRL22_batch;PNRL22_temp];
                    PNWT22_batch = [PNWT22_batch;PNWT22_temp];
                    clear PNRL21_temp
                    clear PNRL22_temp
                    clear PNWT22_temp
                end
                PNRL21_temp_hist = hist(PNRL21_batch,time);
                PNRL22_temp_hist = hist(PNRL22_batch,time);
                PNWT22_temp_hist = hist(PNWT22_batch,time);

                
                PNRL21 = PNRL21+PNRL21_temp_hist;
                PNRL22 = PNRL22+PNRL22_temp_hist;
                PNWT22 = PNWT22+PNWT22_temp_hist;
                sum(PNRL21)
                        
    
                
                if lbatch <buffer_size; %stop reading in batches when we have reached end of .ht3 file.
                    break;
                end                    
            end
            
            PNRL21 = PNRL21(1:bins/1000);
            PNRL22 = PNRL22(1:bins/1000);
            PNWT22 = PNWT22(1:bins/1000);
            time = time(1:bins/1000);

            if ~isprop(obj,file_out_key)
                dummy=addprop(obj,file_out_key)
            end
            obj.(file_out_key)=struct('time',time,'PNRL21',PNRL21,'PNRL22',PNRL22,'PNWT22',PNWT22)
            figure;
            semilogy(time/1000,PNRL21)
            hold on
            semilogy(time/1000,PNRL22)
            semilogy(time/1000,PNWT22)
            

        end
        function obj=get_intensity_correlation(obj,file_key_in,file_key_out,mode,bin_width,time_scale,n_channels);
            %get the intensity correlation function from Toms code.The
            %bin_width is in pulses for t3 data and in ps for t2 data.
            
            %THIS ONLY WORKS FOR MAC OS! WINDOWS OR LINUX USERS HAVE TO
            %FIND A DIFFERENT WAY TO INTERFACE WITH TOMS CODE.
            
            command=['photon_intensity_correlate --file-in ',[obj.pathstr,file_key_in,'.photons_int'], ' --file-out ',[obj.pathstr,file_key_out,'.intensity_corr '],' -m ', mode,' -g 2 -w ', bin_width,' --time-scale ',time_scale,' -c ', num2str(n_channels)]
            system(command)%calls the system command to get intensity_correlation_functions - only on MAC OS
                        
        
        end
        function obj=get_photon_gn_mac(obj,file_key_in,file_key_out,channels,order,time,pulse)
            command=['photon_gn --file-in ',[obj.pathstr,file_key_in,'.photons_int'], ' --file-out ',[obj.pathstr,file_key_out,'.ht3 '],' --mode ', t3,' --order ', order, '--pulse', pulse, ' --time ', time,' --channels ', num2str(channels)]
            

            system(command)%calls the system command to get intensity_correlation_functions - only on MAC OS    

        end
        function obj=get_photon_gn_windows(obj,Data_Path,file_key_in,file_key_out,channels,order,time,pulse)
            % All input values need to be strings, not numbers. Data_Path 
            % is the folder where your data is stored. file key in
            % must be an ht3 file. File key out will determine the folder
            % name of the output file. Channels is the number of detectors
            % used (usually 2 or 4). Order is the correlation order. For
            % two photon events this is 2. Time is actually three numbers
            % separated by commas (no spaces). The first and last should be
            % opposites to create a symmetric function. The determine the
            % time away from the pulse which you extend looking for multi
            % photon events in picoseconds. The center number is the number
            % of time bins which that is divided into and must be an odd
            % number. (ie -400500,801,400500 -- goes +/- 400.5 ns from a 
            % pulse and divides that area into 801 bins). Pulse is similar
            % in structure to time. You need the pulse separation which you
            % extend in either direction. The way the code is written, it
            % looks for the pulse before the number you use, so you must
            % extend an extra 0.5 beyond the pulse you want to go to. To
            % just get center and side peaks you would use '-1.5,3,1.5'
            % which takes three pulses the minus tau, tau = 0 and plus tau
            % peaks. The center number should always be the sum of the
            % absolute value of the plus and minus values.
            command = ['"C:\cygwin64\bin\bash" -c "cd ' Data_Path '; C:/cygwin64/bin/picoquant --file-in ' file_key_in '| C:/cygwin64/bin/photon_gn --mode t3 --channels ' channels ...
                ' --order ' order ' --time ' time ' --pulse ' pulse ' --file-out ' file_key_out '"'];
            system(command)
        end       
        function obj=parse_g2(obj,input_folder,file_key_out,rrate,channels);
            %input_folder is the folder that Thomas's code creates without
            %any extension
            %rrate is repetition rate in Hz
            %Parses the g2 and intensity data structures (originally generated by Tom's
            %correlation code) and converts time units to reasonable values.  In this
            %case, we're parsing the results from a t3 run with time and pulse
            %difference.  However, we're going to assume that the run is
            %pulse-integrated.  That is, each intensity bin is the integral over an
            %entire pulse.  In this case, we clip out the time columns and convert from
            %pulse-units to time units

            %Takes folder from Thomas's code and picks out the g2 file
            file_name_corr = {input_folder};
            file_name_corr = strjoin(file_name_corr,'\');
            file_name_corr = {file_name_corr,'g2.run'};
            file_name_corr = strjoin(file_name_corr,'.');
            file_name_corr = {file_name_corr,'g2'};
            file_name_corr = strjoin(file_name_corr,'\');
            g2 = dlmread(file_name_corr);
            
            %Takes folder from Thomas's code and picks out the intensity
            file_name_corr = {input_folder};
            file_name_corr = strjoin(file_name_corr,'\');
            file_name_corr = {file_name_corr,'g2.run'};
            file_name_corr = strjoin(file_name_corr,'.');
            file_name_corr = {file_name_corr,'intensity'};
            file_name_corr = strjoin(file_name_corr,'\');
            intensity = dlmread(file_name_corr);
            
            %clip out time columns
%             g2 = g2(:,[1:4,7]);

            %timeconversion
            corr_t_units = 1000/rrate; %converts from pulse separation to milliseconds
            intensity_t_units = 1/(rrate); %converts from pulse number to seconds

            %PARSE g2
            %Find out how long the data is.  Each g2 structure contains the 00/01/10/11
            %correlation data.
            l = size(g2,1)/channels^2;

            %Extract the correlation data
            for ii=1:channels^2
                g2_split(:,ii)=g2((ii-1)*l+1:ii*l,5);
                tau_split(:,ii)=0.5*(g2((ii-1)*l+1:ii*l,3)+g2((ii-1)*l+1:ii*l,4))*corr_t_units;
              
            end 
                    
            %Extract a centered-bin time axis for each correlation measurement (they
            %should all be the same, but who knows what really happened, I guess...)
            
            %PARSE intensity
            intensity_t = 0.5*(intensity(:,1)+intensity(:,2))*intensity_t_units;
            
            for ii=1:channels
                intensity_split(:,ii)=intensity(:,ii+2)./((intensity(:,2)-intensity(:,1))*intensity_t_units);             
            end 
            
            intensity_total = sum(intensity_split,2);

            %normalize correlation functions; this will be based on all of the time
            %units being the same.
            totaltime = (intensity(end,1)-intensity(1,1))*corr_t_units; 
            for ii=3:channels+2
                AverageI(ii-2) = sum(intensity(:,ii))/totaltime;
            end
            
            ScalingFactor = prod(AverageI)*(g2(1,4)-g2(1,3))*corr_t_units*totaltime;
            

            count = 0;
            for ii = 1:channels
                for jj = 1:channels
                    count = count+1;
                    g2_norm(:,count) = g2_split(:,count)./(AverageI(ii)*AverageI(jj)*(g2(1:l,4)-g2(1:l,3))*corr_t_units*totaltime);
                end
            end

            %Averaging cross correlations to produce symmetric result
            taucross = tau_split(:,2);
            count = 0;
            g2cross = zeros(l,1);
            for ii = 1:channels
                for jj = 1:channels
                    count = count + 1;
                    if ii == jj
                        continue
                    end
                    g2cross(:) = g2cross(:)+g2_norm(:,count);
                end
            end
            g2cross = g2cross/(channels*(channels-1));


            %plot results for sanity check
            figure
            hold all
            for ii = 1:channels^2
                plot(tau_split(:,ii),g2_split(:,ii))
            end
            title('unnormalized correlation functions')
            
            
            figure
            hold all
            plot(intensity_t,intensity_split)
            title('intensity traces')

            figure
            hold all
            for ii = 1:channels^2
                plot(tau_split(:,ii),g2_norm(:,ii))
            end
            title('normalized correlation functions')
        
           
            %Script that performs the solution-g2 analysis.  It assumes you 
            %have run parse_g2 to populate the workspace.  Furthermore, it requires
            %that the data be a pulse-wise correlation function with symmetric,
            %linearly spaced bins.

            %Parameters
            PCH = 0;
            FCS_Allowoffset=0;
            Plot = 'True';

            %Fit the cross correlation to the FCS governing equations to
            %verify normalization and determine average occupation.
            if FCS_Allowoffset
                OffsetUB = 1.01;
                OffsetLB = .99;
            else
                OffsetUB = 1.00001;
                OffsetLB = .99999;
            end
            %FCS FIT OPTIONS
            fiteqn = 'a1./(1+b1*abs(x))+c1';
            coeffs = {'a1','b1','c1'};
            start = [0.3 0.8 1];
            upper = [10 1000 OffsetUB];
            lower = [0 .001 OffsetLB];
            %weights = 1./taucross((l+3)/2:l);
            weights = 1./taucross((l+11)/2:l);
            fo_ = fitoptions('method','NonlinearLeastSquares','Lower',lower,'Upper',upper,'MaxFunEvals',60000,'MaxIter',4000,'TolFun',1e-008,'TolX',1e-008,'Weights',weights);
            set(fo_,'Startpoint',start);
            ft_ = fittype(fiteqn,'dependent',{'y'},'independent',{'x'},'coefficients',coeffs);  
            %DO FCS FIT
            %FCSfitobject = fit(taucross((l+3)/2:l),g2cross((l+3)/2:l),ft_,fo_);
            FCSfitobject = fit(taucross((l+11)/2:l),g2cross((l+11)/2:l),ft_,fo_);
            FCSfit = FCSfitobject(taucross);
            FCSDiffusionTime = 1/FCSfitobject.b1; %characteristic translational diffusion time in msec
            FCSAvgOccupation = 1/FCSfitobject.a1;
            FCSNormalizationOffset = FCSfitobject.c1-1;
            FCSResidual = mean((g2cross((l+3)/2:l)-FCSfit((l+3)/2:l)).^2);
            %PHENOMENOLOGICAL FIT
            Phenomfitobject = fit(taucross((l+3)/2:(l+3)/2+99),g2cross((l+3)/2:(l+3)/2+99),'poly5');
            Phenomfit = Phenomfitobject(abs(taucross((l+1)/2-100:(l+1)/2+100)));
            PhenomAvgOccupation = 1/(Phenomfit(101)-1);

            if Plot
                %Plot fit
                figure
                subplot(3,3,1:6)
                semilogx(taucross((l+3)/2:l),g2cross((l+3)/2:l));
                hold all
                semilogx(taucross((l+3)/2:l),FCSfit((l+3)/2:l));
                hold off
                xlabel('tau (ms)')
                ylabel('FCS cross-correlation')
                title(strcat('AvgOccupation: ',num2str(FCSAvgOccupation),' Diffusion Time (ms): ',num2str(FCSDiffusionTime),'FCS Offset: ',num2str(FCSNormalizationOffset)))
                subplot(3,3,7:9)
                semilogx(taucross((l+3)/2:l),g2cross((l+3)/2:l)-FCSfit((l+3)/2:l));
                xlabel('tau (ms)')
                ylabel('Residual')

                %plot center results
                figure
                hold all
                plot(taucross((l+1)/2-100:(l+1)/2+100)*1e3,g2cross((l+1)/2-100:(l+1)/2+100));
                plot(taucross((l+1)/2-100:(l+1)/2+100)*1e3,FCSfit((l+1)/2-100:(l+1)/2+100));
                plot(taucross((l+1)/2-100:(l+1)/2+100)*1e3,Phenomfit);
                hold off
                xlabel('tau (us)')
                ylabel('Peak-integrated cross-correlation')
                title(strcat('Phenomenological AvgOccupation: ',num2str(PhenomAvgOccupation)))
            end

            %Calculate quantum yield ratios and stddev under two models (Uses scaling
            %factor variable now added to Parse_g2)
            CenterValue_Phenom = g2cross((l+1)/2);
            CenterStdDev_Phenom = sqrt(CenterValue_Phenom*ScalingFactor)/ScalingFactor;
            SideValue_Phenom = 1+1/PhenomAvgOccupation;
            SideStdDev_Phenom = sqrt(SideValue_Phenom*ScalingFactor)/ScalingFactor;

            CenterValue_FCS = g2cross((l+1)/2);
            CenterStdDev_FCS = sqrt(CenterValue_FCS*ScalingFactor)/ScalingFactor;
            SideValue_FCS = 1+1/FCSAvgOccupation;
            SideStdDev_FCS = sqrt(SideValue_FCS*ScalingFactor)/ScalingFactor;

            PhenomRatio = (CenterValue_Phenom-1)/(SideValue_Phenom-1);
            PhenomStdDev = PhenomRatio*sqrt((CenterStdDev_Phenom/(CenterValue_Phenom-1))^2 + (SideStdDev_Phenom/(SideValue_Phenom-1))^2);
            FCSRatio = (CenterValue_FCS-1)/(SideValue_FCS-1);
            FCSStdDev = FCSRatio*sqrt((CenterStdDev_FCS/(CenterValue_FCS-1))^2 + (SideStdDev_FCS/(SideValue_FCS-1))^2);

            if Plot
            %Report
            disp(strcat('Based on the FCS fit, the quantum yield ratio is: ',num2str(FCSRatio), '.  StdDev is: ',num2str(FCSStdDev)));
            disp(strcat('Based on the phenomenological fit, the quantum yield ratio is: ',num2str(PhenomRatio), '.  StdDev is: ',num2str(PhenomStdDev)));
            end

            if PCH
                %Check for aggregation by examining the intensity profile of the focal
                %volume.
                [nelements,centers] = hist(intensity_01, linspace(min(intensity_01),max(intensity_01),100));

                %Fit with Poissonian model
                fiteqn = 'c1*a1.^(b1*x)*exp(-a1)./gamma(b1*x+1)';
                coeffs = {'a1','b1','c1'};
                start = [10 .001 2000];
                fo_ = fitoptions('method','NonlinearLeastSquares','MaxFunEvals',60000,'MaxIter',4000,'TolFun',1e-008,'TolX',1e-008);
                set(fo_,'Startpoint',start);
                ft_ = fittype(fiteqn,'dependent',{'y'},'independent',{'x'},'coefficients',coeffs);
                %DO FCS FIT
                Intfitobject = fit(centers',nelements',ft_,fo_);
                figure
                hold on
                plot(centers,Intfitobject(centers))
                plot(centers,nelements)
            end
        end
        
        %% C) high level function for analyzing derived properties of the photon stream (fitting FCS traces etc.)
        function obj=fit_auto_corr_FCS_trace(obj,file_key_in,p0,afterpulsing_time)
            
            if nargin<4
                afterpulsing_time=2E6;%given in ps.
            end
            
            %load the auto-correlation function defined by file_key_in
            data=csvread([obj.pathstr,file_key_in,'.intensity_corr']);
            tau=(data(:,3)+data(:,4))./2;
            
            index=[];
            [foo,index(1)]=min((tau-afterpulsing_time).^2); %only treat data after the detector deadtime.
            index(2)=min(find(data(:,5)==0)); %getting index of the first zero value in the correlation function.
            
            tau_select=tau(index(1):index(2)-1);%discarding the last point of tau.
            autocorrelation=data(index(1):index(2)-1,5);%discarding early (afterpulsing) and late (g(2)==0) data. 

            % Define exponential function
            fh = @(x,p) p(1) + p(4).*(1./((1+(x./p(2))).*(1+(p(3).^(-2)).*x./p(2)).^0.5));
            
            % define the error function. this is the function to
            % minimize: you can also use norm or whatever:
            errfh = @(p,x,y) sum((y(:)-fh(x(:),p)).^2);
            
            % an initial guess of the exponential parameters
            if nargin<3
             p0 = [min(autocorrelation) 3E9 1 max(autocorrelation)];
            end
            
            % search for solution
            P = fminsearch(errfh,p0,[],tau_select,autocorrelation);
            residuals=(autocorrelation-fh(tau_select,P));
            fit_curve=fh(tau,P);%extrapolating the curve for all tau.
            tau_limits_for_fit=[min(tau) max(tau)];
            
            %writing the fit results as object properties.
            
            if ~isprop(obj,'FCS_fit')
                dummy=addprop(obj,'FCS_fit');
            end
            
            %wrapping results in a structure. 
            obj.FCS_fit=struct('P',P,'residuals', residuals,'fit_curve',fit_curve,'tau_limits_for_fit',tau_limits_for_fit,'tau_for_fit',tau_select);
           
        end
    end
    end
