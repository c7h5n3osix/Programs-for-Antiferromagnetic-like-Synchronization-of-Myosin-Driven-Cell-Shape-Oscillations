function phase_diagram_E_vs_Eta_Fixed_Final
%% 0. 环境初始化
    % 清理工作区，防止并行计算冲突
    clear; 
    clc;
    
    % 设置全局绘图参数
    set(groot, 'defaultAxesXColor', [0 0 0]);
    set(groot, 'defaultAxesYColor', [0 0 0]);
    set(groot, 'defaultAxesZColor', [0 0 0]);
    set(groot, 'defaultTextColor', [0 0 0]);
    set(groot, 'defaultAxesFontSize', 12);
    
    tic;
    % 启动并行池
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool; 
    end
    
    %% 1. 物理参数定义
    p.k_off   = 1;      
    p.alpha   = 0.01;     
    p.sigma_L = 30;      
    p.K       = 50;     
    
    % --- 固定参数 ---
    const_D     = 1;      % 固定扩散系数 D
    fixed_sigma = 14;     % 【关键】固定活性应力 Sigma，以此绘制 E-Eta 切面
    
    % 导出标量供并行循环使用 (减少通信开销)
    p_Cs = 1 / (1 + p.k_off); 
    p_k_off = p.k_off; 
    p_alpha = p.alpha; 
    p_sigma_L = p.sigma_L; 
    p_K = p.K;
    
    % 绘图限制
    plot_limit_E = 8;    % Y轴 (E) 上限
    
    param_str = sprintf('D=%g, \\sigma_0=%g, \\alpha=%g, k_{off}=%g, K=%g', ...
                        const_D, fixed_sigma, p.alpha, p.k_off, p.K);
    
    %% 2. 扫描参数设置
    resolution = 300; % 分辨率
    
    % X轴: Viscosity (Eta)
    eta_max = 1;
    eta_vec = linspace(0.001, eta_max, resolution)'; 
    
    % 波向量 K 设置 (预计算为行向量)
    k_val_vec = 0.05:0.05:15; 
    k_vec = k_val_vec * pi;   
    K2_vec = (k_vec.^2); 
    
    % 结果预分配
    Boundary_E_RH    = zeros(length(eta_vec), 1); 
    Boundary_E_Delta = zeros(length(eta_vec), 1); 
    
    %% 3. 核心计算循环 (遍历 Eta, 求解临界 E)
    
    limit_val = 10000; 
    
    % 预计算 Ak (Ak 只与 D 有关，不随 Eta 或 E 变化)
    Ak_const = 1 + p_k_off + const_D * K2_vec;
    
    parfor i = 1:length(eta_vec)
        current_eta = eta_vec(i);
        
        % Bk 随 eta 变化
        Bk_vec = 1 + current_eta * K2_vec;
        
        % ==========================================================
        % Part A: Hopf Boundary (解析解法)
        % ==========================================================
        % 寻找最小的 E，使得在 fixed_sigma 下系统稳定。
        % 方程: A(E)*Sigma^2 + B(E)*Sigma + C(E) = 0
        % 转化为 E 的二次方程求解。
        
        term_common = Ak_const + p_K - p_alpha * p_sigma_L;
        
        T = Bk_vec .* term_common;
        U = K2_vec .* term_common;
        V = Ak_const .* Bk_vec * p_K;
        W = p_K * K2_vec .* Ak_const;
        
        m2 = -p_Cs * K2_vec;
        m1 = -p_Cs * p_K * K2_vec;
        
        % 构造 E 的系数
        Qa = m1 .* m2; % Sigma^2 的系数 (不含 E)
        
        Qc2 = U .* K2_vec;
        Qc1 = U .* T + V .* K2_vec - W .* Bk_vec;
        Qc0 = V .* T;
        
        Qb1 = m1 .* K2_vec + m2 .* U;
        Qb0 = m1 .* T + m2 .* V;
        
        Sig = fixed_sigma;
        
        % E 的二次方程: (Qc2)E^2 + (Qb1*Sig + Qc1)E + (Qa*Sig^2 + Qb0*Sig + Qc0) = 0
        Poly_A = Qc2;
        Poly_B = Qb1 * Sig + Qc1;
        Poly_C = Qa * Sig^2 + Qb0 * Sig + Qc0;
        
        Discriminant_E = Poly_B.^2 - 4 * Poly_A .* Poly_C;
        valid_k = Discriminant_E >= 0;
        
        max_required_E = 0; 
        
        if any(valid_k)
            E_root1 = (-Poly_B(valid_k) + sqrt(Discriminant_E(valid_k))) ./ (2 * Poly_A(valid_k));
            E_root2 = (-Poly_B(valid_k) - sqrt(Discriminant_E(valid_k))) ./ (2 * Poly_A(valid_k));
            
            all_roots = [E_root1, E_root2];
            real_roots = all_roots(all_roots > 0);
            
            if ~isempty(real_roots)
                % 取所有不稳定 k 中最大的那个 E 作为边界
                max_required_E = max(real_roots);
            end
        end
        Boundary_E_RH(i) = max_required_E;
        
        
        % ==========================================================
        % Part B: Delta Boundary (二分法求解 E)
        % ==========================================================
        % 目标：找到 E，使得系统的 Collapse 临界活性值等于 fixed_sigma
        
        low = 0; 
        high = plot_limit_E * 2; % 搜索范围
        
        % 计算端点处的 Critical Sigma
        val_low  = get_critical_sigma(low,  K2_vec, Ak_const, Bk_vec, p_Cs, p_K, p_alpha, p_sigma_L);
        val_high = get_critical_sigma(high, K2_vec, Ak_const, Bk_vec, p_Cs, p_K, p_alpha, p_sigma_L);
        
        % 目标是 val(E) - fixed_sigma = 0
        % 通常 E 越大，Critical Sigma 越大 (越难 Collapse)
        
        found_root = limit_val;
        
        if val_low > fixed_sigma
            % 即使 E=0，临界值也比当前活性高 -> 永远不会 Collapse
            found_root = 0;
        elseif val_high < fixed_sigma
            % 即使 E 很大，临界值也比当前活性低 -> 总是 Collapse (或者超出绘图范围)
            found_root = limit_val;
        else
            % 二分搜索
            for iter = 1:25
                mid = (low + high) / 2;
                val_mid = get_critical_sigma(mid, K2_vec, Ak_const, Bk_vec, p_Cs, p_K, p_alpha, p_sigma_L);
                
                if val_mid > fixed_sigma
                    % 中点的临界值太高了，说明 E 给大了 (太稳了)
                    high = mid;
                else
                    % 中点的临界值太低了，说明 E 不够，容易 Collapse
                    low = mid;
                end
            end
            found_root = high;
        end
        Boundary_E_Delta(i) = found_root;
        
    end
    
    %% 4. 绘图与区域填充
    
    % 数据截断
    Boundary_E_RH(Boundary_E_RH > plot_limit_E) = plot_limit_E;
    Boundary_E_Delta(Boundary_E_Delta > plot_limit_E) = plot_limit_E;
    
    % 定义区域边界
    % Stable: 需要 E 足够大，大于 RH 要求且大于 Delta 要求
    Upper_Boundary = max(Boundary_E_RH, Boundary_E_Delta);
    
    figure('Color', 'white', 'Position', [100, 100, 750, 600]);
    hold on; box off;
    
    Color_Stable   = [0.92, 0.96, 1.00]; 
    Color_Osc      = [0.60, 0.80, 0.95]; 
    Color_Collapse = [0.25, 0.40, 0.60]; 
    
    % 1. 填充 Stable (上方区域)
    fill([eta_vec; flipud(eta_vec)], [Upper_Boundary; ones(size(eta_vec))*plot_limit_E], Color_Stable, ...
        'EdgeColor', 'none', 'DisplayName', 'Stable');
    
    % 2. 填充 Collapse (下方区域，E < Delta Boundary)
    fill([eta_vec; flipud(eta_vec)], [zeros(size(eta_vec)); flipud(Boundary_E_Delta)], Color_Collapse, ...
        'EdgeColor', 'none', 'DisplayName', 'Collapse');
        
    % 3. 填充 Oscillatory (中间夹层: Delta < E < RH)
    % 只有当 RH 边界 高于 Delta 边界时，才存在振荡区
    Band_Bottom = Boundary_E_Delta;
    Band_Top = max(Boundary_E_Delta, Boundary_E_RH);
    
    fill([eta_vec; flipud(eta_vec)], [Band_Bottom; flipud(Band_Top)], Color_Osc, ...
        'EdgeColor', 'none', 'DisplayName', 'Oscillatory');
    
    % 绘制线条
    plot(eta_vec, Boundary_E_RH, 'b--', 'LineWidth', 2, 'DisplayName', 'Hopf Boundary');
    plot(eta_vec, Boundary_E_Delta, 'r-', 'LineWidth', 2, 'DisplayName', '\Delta=0 Boundary');
    
    % 装饰
    xlim([0, eta_max]); 
    ylim([0, plot_limit_E]);
    xlabel('Viscosity \eta', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Elastic Modulus E', 'FontSize', 14, 'FontWeight', 'bold');
    
    title({'Phase Diagram (E vs \eta)', ['\fontsize{10}\color[rgb]{0.3,0.3,0.3}' param_str]}, ...
          'Interpreter', 'tex', 'FontWeight', 'bold');
      
    legend('Location', 'best');
    set(gca, 'Layer', 'top', 'LineWidth', 1.5, 'FontSize', 12, 'TickDir', 'in');
    
    % 标注
    text(eta_max*0.5, plot_limit_E*0.9, 'Stable', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0 0.2 0.5], 'HorizontalAlignment', 'center');
    text(eta_max*0.5, 1, 'Collapse', 'FontSize', 16, 'FontWeight', 'bold', 'Color', 'white', 'HorizontalAlignment', 'center');
    
    toc;
end

%% 辅助函数: 计算给定 E 下的临界 Collapse Sigma
function val = get_critical_sigma(E_val, K2, Ak, Bk, Cs, pK, alpha, sigL)
    % 该函数计算在特定 E 和 eta (隐含在 Bk 中) 下，
    % 系统发生实根分岔(Collapse)的最小 Sigma 值。
    
    m2 = -Cs * K2; 
    
    term_common = Ak + pK - alpha * sigL;
    c2 = Bk .* term_common + E_val * K2;
    
    m1 = -Cs * pK * K2;
    c1 = E_val * K2 .* term_common + Ak .* Bk * pK;
    
    c0 = E_val * pK * K2 .* Ak;
    c3 = Bk;
    
    min_sigma = 10000; % 初始大值
    
    % 遍历 K (由于 roots 无法向量化，这里使用循环)
    for j = 1:length(K2)
        mm2 = m2(j); cc2 = c2(j);
        mm1 = m1(j); cc1 = c1(j);
        cc0 = c0(j); cc3 = c3(j);
        
        % 手动展开卷积以提高计算速度
        
        % T3 (x^4系数)
        % pA2 = [mm2, cc2], pA2^2 = [mm2^2, 2mm2cc2, cc2^2]
        A2_sq_1 = mm2^2; A2_sq_2 = 2*mm2*cc2; A2_sq_3 = cc2^2;
        A1_sq_1 = mm1^2; A1_sq_2 = 2*mm1*cc1; A1_sq_3 = cc1^2;
        
        % conv(A2^2, A1^2)
        T3_1 = A2_sq_1 * A1_sq_1;
        T3_2 = A2_sq_1 * A1_sq_2 + A2_sq_2 * A1_sq_1;
        T3_3 = A2_sq_1 * A1_sq_3 + A2_sq_2 * A1_sq_2 + A2_sq_3 * A1_sq_1;
        T3_4 = A2_sq_2 * A1_sq_3 + A2_sq_3 * A1_sq_2;
        T3_5 = A2_sq_3 * A1_sq_3;
        
        % T2 (x^3系数)
        termT2 = 4*cc0;
        T2_1 = termT2 * (mm2^3);
        T2_2 = termT2 * (3*mm2^2*cc2);
        T2_3 = termT2 * (3*mm2*cc2^2);
        T2_4 = termT2 * (cc2^3);
        
        % T4 (x^3系数)
        termT4 = 4*cc3;
        T4_1 = termT4 * (mm1^3);
        T4_2 = termT4 * (3*mm1^2*cc1);
        T4_3 = termT4 * (3*mm1*cc1^2);
        T4_4 = termT4 * (cc1^3);
        
        % T1 (x^2系数)
        termT1 = 18 * cc3 * cc0;
        T1_1 = termT1 * (mm2*mm1);
        T1_2 = termT1 * (mm2*cc1 + cc2*mm1);
        T1_3 = termT1 * (cc2*cc1);
        
        % T5 (x^0系数)
        T5 = 27 * cc3^2 * cc0^2;
        
        % 合并系数
        C1 = T3_1;
        C2 = T3_2 - T2_1 - T4_1;
        C3 = T3_3 - T2_2 - T4_2 + T1_1;
        C4 = T3_4 - T2_3 - T4_3 + T1_2;
        C5 = T3_5 - T2_4 - T4_4 + T1_3 - T5;
        
        rr = roots([C1 C2 C3 C4 C5]);
        valid_r = rr(imag(rr)==0 & real(rr)>0);
        
        if ~isempty(valid_r)
            % 对于当前 k，该 E 值导致 Collapse 的临界 Sigma 是 max(roots)
            k_crit = max(valid_r);
            if k_crit < min_sigma
                min_sigma = k_crit;
            end
        end
    end
    val = min_sigma;
end