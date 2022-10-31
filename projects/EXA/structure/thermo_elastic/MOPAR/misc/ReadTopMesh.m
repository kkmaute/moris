function [fecoord,fetopo,neib,ebtyp]=ReadTopMesh(topfile)

maxnodes=500000;
maxelem=1000000;
maxblk=10;

fid=fopen(topfile,'r');

key=0;
numnodes=0;
numele=0;
numblk=0;

nodewarn=0;
elemwarn=0;

fecoord = zeros(maxnodes,3);
fetopo  = zeros(maxelem,8);
neib    = zeros(maxblk,1);
ebtyp   = cell(maxblk,1);

while 1
    
    % read line
    line=fgetl(fid);
    
    % check of EOF
    if line<0;break;end
    
    % check for key word: element or node
    if strncmp('Nodes',line,5) > 0; 
        key=10;
    end
    if strncmp('Elements',line,8) > 0; 
        key=20; 
        numblk=numblk+1;
    end
    
    switch key
        
        case 1  % nodal coordinates
            
            nddata=sscanf(line,'%d %e %e %e');
            fecoord(nddata(1),:)=nddata(2:4);
            numnodes=max(numnodes,nddata(1));
            
            if numnodes > maxnodes && nodewarn==0
                fprintf('warning numnodes exceeds maxnodes\n');
                nodewarn=1;
            end
            
        case 2 % elements
        
            neib(numblk)=neib(numblk)+1;
            eldata=sscanf(line,'%d %d');

            numele=max(numele,eldata(1));

            if numelem > maxelem && elemwarn==0
                fprintf('warning numnodes exceeds maxnodes\n');
                elemwarn=1;
            end
            
            switch eldata(2)
                
                case 4 % 3-node triangle
                   ebtyp{numblk}='tri3';
                   eldata=sscanf(line,'%d %d %d %d %d');
                   fetopo(eldata(1),1:4)=eldata(2:5);
                   
                case 5 % 4-node tet
                   ebtyp{numblk}='tet4';
                   eldata=sscanf(line,'%d %d %d %d %d %d');
                   fetopo(eldata(1),1:5)=eldata(2:6);
                   
                otherwise
                    error('incorrect element typ');
            end
            
    end
    
    if key>9; key=key/10; end
end

fclose(fid);

% trim arrays
fecoord = fecoord(1:numnodes,:);
fetopo  = fetopo(1:numele,:);
neib    = neib(1:numblk,:);
ebtyp   = ebtyp(1:numblk,:);
