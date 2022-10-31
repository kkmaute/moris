deps=1e-10;

for ie=1:1600
    load(['rjele_' num2str(ie)],'re','je','me');
    re_new=re;
    je_new=je;
    me_new=me;

    load(['/home/maute/Work/TransFem/Studies/Bender_XFEM/rjele_' num2str(ie)],'re','je','me');
    re_old=re;
    je_old=je;
    me_old=me;
    
    rediff=norm(re_new-re_old)/max(norm(re_new+re_old),eps);
    jediff=norm(je_new-je_old)/max(norm(je_new+je_old),eps);
    mediff=norm(me_new-me_old)/max(norm(me_new+me_old),eps);

%     rediff=max(max(abs(re_new-re_old)))/max(norm(re_new+re_old),eps);
%     jediff=max(max(abs(je_new-je_old)))/max(norm(je_new+je_old),eps);
%     mediff=max(max(abs(me_new-me_old)))/max(norm(me_new+me_old),eps);

    if rediff>deps || jediff>deps || mediff>deps
        fprintf('ie=%d  re=%e  je=%e  me=%e\n',ie,rediff,jediff,mediff);
    end
    
end
