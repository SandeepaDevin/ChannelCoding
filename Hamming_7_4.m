EbNo=-4:2:8;
for i=1:length(EbNo)
    EbNo1=EbNo(i);BER=0;errors=0;data=0;
    disp(EbNo(i));
    Nodb=10.^(EbNo1/10);
while errors<100
    genbit=randi([0 1],1,11);
    %% Hamming encoding process
     windowq=length(genbit)/11;
     fbitwind=reshape(genbit,11,windowq)';
     P= [1 1 0 0 ;0 1 1 0 ;0 0 1 1 ;1 0 1 0 ;1 0 0 1 ;0 1 0 1 ;1 1 1 0 ;0 1 1 1 ;1 0 1 1 ;1 1 0 1 ;1 1 1 1; ]; % Coeffcient Matrix
     generator = [P eye(11)];
     out=mod(fbitwind*generator,2);
     out=reshape(out',numel(out),1)';
     %% bpsk mod and transmission
     bpsk=pskmod(out,2);
     noice=normrnd(0,0.1^(0.5));
     Eb=(0.2)*10^(EbNo1/10);
     trans=(sqrt(Eb)*bpsk)+noice;
     demod=pskdemod(trans,2);
     %% hamming decoding process
     NoofWin=numel(demod)/15;
     Hammarray = zeros(round(length(genbit)*(15/11)), 1);
     sbwin=reshape(demod,15,NoofWin)';
     h=[eye(4) P' ];
     Syn= mod(sbwin*h',2); syntbl=syndtable(h);
     indexes=Syn*[8;4;2;1];
     err=zeros(NoofWin,15);
     for g = 1:NoofWin
     err(g,:) = syntbl(indexes(g)+1,:);
     end
     Hammarray = mod(sbwin + err, 2); 
     Hammarray(:,1:4) = [];                            
     outpuut= (reshape(Hammarray',length(genbit), 1))';
     %% error detection
    for j = 1:length(genbit)
     if outpuut(j)~=genbit(j)
         errors = errors+1;        
     end
    end
    data=data+11;
end
 BER=errors/data;
 BER_th=qfunc(sqrt(2*(15/11)*Nodb));
 disp(BER)
 disp(BER_th)
end    
