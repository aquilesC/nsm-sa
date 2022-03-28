
%% scatter plot
for i=1:length(collection.iOC_mean)
                scatter(collection.iOC_mean(i),collection.Deff_mean(i),'MarkerFaceAlpha',collection.N(i)/max(collection.N),...
                    'MarkerFaceColor','black',...
                    'MarkerEdgeColor','none'); hold on
end