function scatterplot_with_correlation(data,x_vars,y_vars, sig_level,savname,varargin)
% sig_level should be a scalar between 0 and 1 representing the desired significance level for color coding
% scatterplot_with_correlation(tdata,col,col, 0.05,savname,...
%    20,'MarkerFaceColor',LC(2,:),'MarkerEdgeColor','none')

m=length(y_vars);
n=length(x_vars);

% Calculate aspect ratio
aspect_ratio = n / m;

% Get display size
screen_size = get(groot,'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);

% Calculate figure size based on aspect ratio and display size
fig_width = 0.5 * screen_width; % example value
fig_height = fig_width / aspect_ratio/1.3 ;

% Set figure position and size
fig_xpos = 0.1 * screen_width; % example value
fig_ypos = 0.1 * screen_height; % example value
fig_position = [fig_xpos, fig_ypos]; % calculated value
fig_size = [fig_width, fig_height]; % calculated value

fig = figure;
fig.Position = [fig_position, fig_size];
set(gcf,'Color','w')

tiledlayout(m,n,TileSpacing="tight",Padding="tight");
ax=zeros(m,n);
for i =1:m
    for j =1:n
        ax(i,j)=nexttile;
    end
end

for i =1:m
    for j =1:n
        x=data.(x_vars{j});y= data.(y_vars{i});
        if strcmpi(x_vars{j},y_vars{i})
            histogram(ax(i,j),x,Normalization="pdf",EdgeColor='none')
        else
            scatter(ax(i,j),x,y,varargin{:});hold(ax(i,j),'on')
            [corr_coeff,p_value] = corr(x,y);
            if p_value < sig_level
                color = 'r';
            else
                color = 'k';
            end
            xx = min(x):0.01:max(x);
            lm=fitlm(x,y);
            yy = predict(lm,xx');
            plot(ax(i,j),xx, yy, 'k--','LineWidth',1.5)
            p1=plot(ax(i,j),nan,nan,'Linestyle', 'none', 'Marker', 'none');
            lh=legend(ax(i,j),p1,sprintf('%.2f', corr_coeff));
            lh.Location='best';lh.Box='off';lh.TextColor=color;
            lh.FontSize=8;
        end
    end
end

for i = 1:m
    ylabel(ax(i,1),y_vars(i), 'Interpreter','latex')
end

for  j =1:n
    xlabel(ax(m,j),x_vars(j), 'Interpreter','latex')
end

if m>=1 && n>=2
    set(ax(1:m-1,:),'XTickLabel',[])
    set(ax(:,2:n),'YTickLabel',[])
end


set(ax, 'Fontsize', 9, 'LineWidth', 0.45,'box','on', ...
    'Xcolor', [1, 1, 1]*0.25, 'Ycolor', [1, 1, 1]*0.25,...
    'Xgrid','off','Ygrid','off')

% export_fig(fig,"results/"+savname,'-r300')

if ~isempty(savname)
    exportgraphics(fig,"results/"+savname,Resolution=300)
end


end
