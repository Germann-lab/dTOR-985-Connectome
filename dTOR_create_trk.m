% Script to generate a .trk file containing MNI space streamlines, as specified by an input .mat file containing streamline weight data
% Only streamlines with weights above a certain threshold will be generated
% Suresh Joel, General Electric Global Research, May 2018

% Load the tract data
baseFolder='/PATH_TO_MAT_FILE'; % Provide path for input .mat file
load(fullfile(baseFolder,'INPUT_MAT_FILE.mat')); % Provide filename for input .mat file

% Load the full dTOR connectome file
inFile = '/PATH_TO_dTOR_CONNECTOME_FILE/dTOR_985.mat'; % Insert path to dTOR connectome file
load(inFile);

nfibers=max(fibers(:,4));

labs = fibers(:,4); % Fiber numbers/labels
dlabs = diff(labs); % Mark indexes where the labels change (assuming increasing order of labels)
ch_ind_labs = find(dlabs); % Find the marked indices
ch_ind_labs = [ch_ind_labs; length(fibers(:,4))];
fibers_cell={};
fib_ind_start=1;
for iFib=1:nfibers
   % tic
    a=fibers(fib_ind_start:ch_ind_labs(iFib),1:3);
    fibers_cell{iFib}=a;
  %  disp(['Fiber ',num2str(iFib),' done']);
    fib_ind_start=ch_ind_labs(iFib)+1;
   % toc
end

%save(outFile,'fibers_cell','-v7.3');

% Specify the filename for the output file
outfile = fullfile(baseFolder,'OUTPUT_FILENAME.trk'); % Provide filename for input .trk file

% Set the streamline weight threshold for the output file (only streamlines with weights larger than this value will be retained)
wt_thresh=0; % Replace '0' with the threshold value to be used, as appropriate
h_th=max(fibers_wt);
l_th=wt_thresh;

g_fib_id=find(fibers_wt>wt_thresh);
nfib=length(g_fib_id);
disp(nfib)

% Create .trk file header
id=['TRACK',char(0)]; %id (6 char)
dim=[181,217,181]; %dim (short)
sz=[1.0,1.0,1.0]; %size (float)
org=[ 0 0 0]; %origin (float)

n_s=3; %n_scalars (short)
n_s_names=['Red',char(zeros(1,17)), ...
    'Green',char(zeros(1,15)), ...
    'Blue',char(zeros(1,16)), ...
    char(zeros(1,140))]'; %scalar name (200 char)
n_p=1; %n_properties (short)
n_p_names=['DBS_Weight',char(zeros(1,190))]'; %property name (200 char)
%vox_to_ras=[1 0 0 90; 0 1 0 126;  0 0 1 72; 0 0 0 1]; %reshape(eye(4),[],1); %vox_to_ras (float)
%origin = [90,126,72];
vox_to_ras=[0 0 0 0; 0 0  0 0;  0 0 0 0; 0 0 0 0]; %reshape(eye(4),[],1); %vox_to_ras (float)
res=char(zeros(444,1)); %reserved (444 char)
vox_ord=['LAS',char(0)]; %voxel order (4 char)
pad2=['LPS',char(0)];%char(zeros(4,1)); %pad2 (4 char)
img_orn_pt=[1 0 0 0 1 0]; %img_orn_pt (6 float)
pad1=char(zeros(2,1)); %pad1 (2 char)
invrt=char([0 0 0]); %invert x y z (3 uchar)
swp=char([0 0 0]); %swap xy yz zx (3 uchar)
n_tracks=nfib; %n_tracks (int32)
ver=2; %version (int32)
hdr_sz=1000; %header size (int32)

% Open the .trk file
fid=fopen(outfile,'W');

% Write header
fwrite(fid,id,'char');
fwrite(fid,dim,'int16');
fwrite(fid,sz,'float');
fwrite(fid,org,'float');
fwrite(fid,n_s,'int16');
fwrite(fid,n_s_names,'char');
fwrite(fid,n_p,'int16');
fwrite(fid,n_p_names,'char');
fwrite(fid,reshape(vox_to_ras',[],1),'float');
fwrite(fid,res,'char');
fwrite(fid,vox_ord,'char');
fwrite(fid,pad2,'char');
fwrite(fid,img_orn_pt,'float');
fwrite(fid,pad1,'char');
fwrite(fid,invrt,'uchar');
fwrite(fid,swp,'uchar');
fwrite(fid,n_tracks,'int32');
fwrite(fid,ver,'int32');
fwrite(fid,hdr_sz,'int32');


% Set colour scheme for .trk file and write data to file
aut=summer; % Specify the colour scheme with which to display streamline weights (refer to https://www.mathworks.com/help/matlab/ref/colormap.html)
for i=1:length(g_fib_id)
    %fib_pts=fibers(fibers(:,4)==g_fib_id(i),1:3);
    fib_pts = fibers_cell{g_fib_id(i)};
    n=size(fib_pts,1);
    fwrite(fid,n,'int32');
    for j=1:n
        %dimxy(3)=2*(fib_pts(j,3)+(origin(3).*orns(iOrn,3)));
        %fwrite(fid,fib_pts(j,:).*[-1,-1,1]+(origin.*orns(iOrn,:)),'float');
        fwrite(fid,fib_pts(j,:).*[-1,-1,1],'float');
        col=round(aut(round((63/(h_th-l_th))*((fibers_wt(g_fib_id(i)))-l_th))+1,:)*255);
        fwrite(fid,col,'float');
    end;
    fwrite(fid,fibers_wt(g_fib_id(i)),'float');
    clear fib_pts n;
end;

% Close the file
fclose(fid);