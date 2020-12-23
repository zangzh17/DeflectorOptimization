classdef Material
    properties
        n_coeffs {mustBeNumeric} = []
        n_model = 'sellm'
        n_order = 4
        k_coeffs {mustBeNumeric} = []
        k_model = 'poly'
        k_order = 6
        name
    end
    methods
        function obj = Material(name,varargin)
            p = inputParser;
            addRequired(p,'name');
            addOptional(p,'n_coeffs',[]);
            parse(p,name,varargin{:});
            obj.name = p.Results.name;
            obj.n_coeffs = p.Results.n_coeffs;
        end
        % 重载索引函数
        function varargout = subsref(obj,S)
            % 字符串索引且是()类型的索引，则重载
            if strcmp(S(1).type, '()')
                if ischar(S(1).subs{1})
                    name_list = cell(size(obj));
                    [name_list{:}] = obj(:).name;
                    if length(S)>1
                        [varargout{1:nargout}] = builtin('subsref',obj((cellfun(@(x) strcmp(x,S(1).subs{1}),name_list))), S(2:end));
                    else
                        [varargout{1:nargout}] = obj((cellfun(@(x) strcmp(x,S.subs{1}),name_list)));
                    end
                    return
                else
                    % 其他情况不重载，调用内置函数
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                end
            else
                % 其他情况不重载，调用内置函数
                [varargout{1:nargout}] = builtin('subsref',obj,S);
            end
        end
        
        function n = getN(obj,lambda)
            % 输入折射率对象和波长（um），输出折射率数据
            switch obj.n_model
                case 'sellm'
                    n = real(Material.sellm(obj.n_coeffs,lambda));
                case 'poly'
                    n = real(polyval(obj.n_coeffs,lambda));
            end
        end
        
        function n = getIndex(obj,lambda)
            % 输入折射率对象和波长（um），输出复折射率数据
            switch obj.n_model
                case 'sellm'
                    n = Material.sellm(obj.n_coeffs,lambda);
                case 'poly'
                    n = polyval(obj.n_coeffs,lambda);
            end
            if ~isempty(obj.k_coeffs)
                switch obj.k_model
                    case 'sellm'
                        k = Material.sellm(obj.k_coeffs,lambda);
                    case 'poly'
                        k = polyval(obj.k_coeffs,lambda);
                end
                n = n +1j*k;
            end
        end
        
        function obj = dataFit(obj,lambda_n,varargin)
            % 调用 inputParser 以创建解析器对象
            p = inputParser;
            addRequired(p,'obj');
            addRequired(p,'lambda_n');
            addOptional(p,'lambda_k',[]);
            addParameter(p,'n_model','sellm');
            addParameter(p,'n_order',4);
            addParameter(p,'k_model','poly');
            addParameter(p,'k_order',6);
            addParameter(p,'w_range',[0,inf]);
            addParameter(p,'fig',[]);
            parse(p,obj,lambda_n,varargin{:});
            % 参数拟合
            % 拟合模型选择
            obj.n_model = p.Results.n_model;
            obj.k_model = p.Results.k_model;
            obj.n_order = p.Results.n_order;
            obj.k_order = p.Results.k_order;
            
            % 折射率拟合
            ind = (~isnan(p.Results.lambda_n(:,1))) &...
                p.Results.lambda_n(:,1)>p.Results.w_range(1) &...
                p.Results.lambda_n(:,1)<p.Results.w_range(2);
            lambda_n = p.Results.lambda_n(ind,:);
            switch obj.n_model
                case 'sellm'
                    % 拟合算法选择，可以选其中一个method{1}~method{4}
                    method = {'andrews','cauchy','bisquare','logistic',[]};
                    opts = statset('nlinfit');
                    opts.RobustWgtFun = method{obj.n_order};
                    coeff0 = [1 1.8 0.1 1 0.12]; % 优化初始值
                    c = nlinfit(lambda_n(:,1), lambda_n(:,2), @Material.sellm , coeff0, opts);
                case 'poly'
                    % 拟合阶数
                    c = polyfit(lambda_n(:,1), lambda_n(:,2), obj.n_order);
            end
            obj.n_coeffs = c;
            % 折射率虚部拟合
            if ~isempty(p.Results.lambda_k)
                ind = (~isnan(p.Results.lambda_k(:,1))) &...
                    p.Results.lambda_k(:,1)>p.Results.w_range(1) &...
                    p.Results.lambda_k(:,1)<p.Results.w_range(2);
                lambda_k = p.Results.lambda_k(ind,:);
                switch obj.k_model
                    case 'sellm'
                        % 拟合算法选择，可以选其中一个method{1}~method{4}
                        method = {'andrews','cauchy','bisquare','logistic',[]};
                        opts = statset('nlinfit');
                        opts.RobustWgtFun = method{obj.k_order};
                        coeff0 = [1 1.8 0.1 1 0.12]; % 优化初始值
                        c = nlinfit(lambda_k(:,1), lambda_k(:,2), @Material.sellm , coeff0, opts);
                    case 'poly'
                        % 拟合阶数
                        c = polyfit(lambda_k(:,1), lambda_k(:,2), obj.k_order);
                end
                obj.k_coeffs = c;
            else
                obj.k_coeffs = [];
            end
            % 绘图
            if ~isempty(p.Results.fig)
                figure(p.Results.fig)
            end
            lambda = linspace(min(lambda_n(:,1)),max(lambda_n(:,1)),100);
            t = tiledlayout('flow','TileSpacing','compact');
            title(t,obj.name)
            nexttile
            plot(lambda_n(:,1),lambda_n(:,2),'.')
            hold on
            plt = plot(lambda,real(getIndex(obj,lambda)),'-','LineWidth',0.5);
            legend(plt, ['model:' obj.n_model ' ord: ' num2str(obj.n_order)])
            ylabel('n')
            if ~isempty(obj.k_coeffs)
                nexttile
                plot(lambda_k(:,1),lambda_k(:,2),'.')
                hold on
                plt = plot(lambda,imag(getIndex(obj,lambda)),'-','LineWidth',0.5);
                legend(plt, ['model:' obj.k_model ' ord: ' num2str(obj.k_order)])
                ylabel('k')
            end
        end
        function showN(obj,fig,w_range)
            if nargin==1
                figure
                w_range = [0.5,2];
            elseif nargin==2
                figure(fig);
                w_range = [0.5,2];
            else
                figure(fig);
            end
            lambda = linspace(min(w_range),max(w_range),100);
            plot(lambda,real(getIndex(obj,lambda)),'-','LineWidth',1.5);
            ylabel('n')
            title(obj.name)
        end
        function showIndex(obj,fig,w_range)
            if nargin==1
                figure
                w_range = [0.5,2];
            elseif nargin==2
                figure(fig);
                w_range = [0.5,2];
            else
                figure(fig);
            end
            lambda = linspace(min(w_range),max(w_range),100);
            if ~isempty(obj.k_coeffs)
                yyaxis left
            end
            plot(lambda,real(getIndex(obj,lambda)),'-','LineWidth',1.5);
            ylabel('n')
            if ~isempty(obj.k_coeffs)
                yyaxis right
                plot(lambda,imag(getIndex(obj,lambda)),'-','LineWidth',1.5);
                ylabel('k')
            end
            title(obj.name)
        end
    end
    methods (Static)
        function nn = sellm(coeff,lambda)
            % sellm模型定义
            % 输入sellm模型参数c，输出折射率数据
            %sel Calculate sellmier eqs
            %   n(l)^2 = A + B1*(x^2)/(x^2-C1^2) + B2*(x^2)/(x^2-C2^2) + ...,
            % where x is the wavelength in um
            L=lambda.^2;
            C = coeff(3:2:end);%odd numbers
            B = coeff(2:2:end); % even numbers
            ssum = 0;
            for ii=1:length(C)
                ssum = ssum + B(ii).*L./(L-C(ii));
            end
            nn=(sqrt(coeff(1) + ssum));
        end
    end
end