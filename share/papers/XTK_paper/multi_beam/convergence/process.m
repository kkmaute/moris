%% Script to check rank of the linear system
close all;
clear;
clc;

%% Read values from .csvs

NoEnrich = readmatrix('Enrichment_false/errors.csv');
WithEnrich = readmatrix('Enrichment_true/errors.csv');
MatSize = size( NoEnrich );
nVals = MatSize(1);

NoEnrichL2Ref = NoEnrich(2);
NoEnrichH1sRef = NoEnrich(3);
WithEnrichL2Ref = WithEnrich(2);
WithEnrichH1sRef = WithEnrich(3);

nDofNoEnrich = NoEnrich(2:nVals,1);
nDofWithEnrich = WithEnrich(2:nVals,1);

RelL2NoEnrich = sqrt( NoEnrich(2:nVals,2)./NoEnrich(1,2) );
RelH1sNoEnrich = sqrt( NoEnrich(2:nVals,3)./NoEnrich(1,3) );

RelL2WithEnrich = sqrt( WithEnrich(2:nVals,2)./WithEnrich(1,2) );
RelH1sWithEnrich = sqrt( WithEnrich(2:nVals,3)./WithEnrich(1,3) );

%% Plot the convergence

% plot L2-convergence
figure;
loglog(nDofNoEnrich,RelL2NoEnrich,'-o');
hold on;
loglog(nDofWithEnrich,RelL2WithEnrich,'-o');
grid on;
legend('No Enrich.','With Enrich.');
xlabel('n_{DOF}');
ylabel('Rel. L2');

% plot H1s-convergence
figure;
loglog(nDofNoEnrich,RelH1sNoEnrich,'-o');
hold on;
loglog(nDofWithEnrich,RelH1sWithEnrich,'-o');
grid on;
legend('No Enrich.','With Enrich.');
xlabel('n_{DOF}');
ylabel('Rel. H1s');

%% Write to CSV for later use

NoEnrichWriteMat = zeros(nVals:3);
NoEnrichWriteMat(:,1) = nDofNoEnrich;
NoEnrichWriteMat(:,2) = RelL2NoEnrich;
NoEnrichWriteMat(:,3) = RelH1sNoEnrich;

WithEnrichWriteMat = zeros(nVals:3);
WithEnrichWriteMat(:,1) = nDofWithEnrich;
WithEnrichWriteMat(:,2) = RelL2WithEnrich;
WithEnrichWriteMat(:,3) = RelH1sWithEnrich;

writematrix( NoEnrichWriteMat, 'errors_non_enriched.csv' );
writematrix( WithEnrichWriteMat, 'errors_enriched.csv' );

