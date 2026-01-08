#this is for V = 1.5
import matplotlib

matplotlib.use('TkAgg')
import numpy as np
import time
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing as mp
from collections import deque  # <--- 1. 引入 deque
# from numba import jit # <--- 1. 导入 jit
#
# @jit(nopython=True)
def AdmissionControl(lam, q, p_R, p_L, nR, nC, nL, V, gamma, beta):
    # Simulation parameters
    T = 500
    R_w = 0.2
    T_w = T * R_w
    # gamma = 0.7
    # beta = 0.5
    mu = 1
    # V = 2
    c = 0.5
    lmd_R = lam * gamma * q
    lmd_L = lam * (1 - gamma) * p_L
    # init
    t = 0.0
    ta_R = np.random.exponential(1.0 / lmd_R) if lmd_R > 0 else np.inf
    ta_L = np.random.exponential(1.0 / lmd_L) if lmd_L > 0 else np.inf
    td = np.inf
    TravelTimes = []
    TravelQueueInfo = []
    OrderQueueInfo = []

    Na_R = Na_L = Nd = Nc = 0
    cumTH = 0
    n_o = 0
    n_t = 0

    Utility_remote = []
    Utility_local = []
    ArrTime_L = []
    ArrTime_R = []

    while t < T:
        # determine next event robustly
        next_travel = min(TravelTimes) if TravelTimes else np.inf
        events = [ta_R, ta_L, td, next_travel]
        #events = np.array([ta_R, ta_L, td, next_travel])
        idx = int(np.argmin(events))
        next_event_time = events[idx]

        if next_event_time == np.inf:
            break

        t = next_event_time

        # 0: remote arrival
        if idx == 0:
            ta_R = t + np.random.exponential(1.0 / lmd_R) if lmd_R > 0 else np.inf
            if n_o < nR:
                Na_R += 1
                ArrTime_R.append(t)
                n_o += 1
                OrderQueueInfo.append(Na_R)  # positive id for remote
                if n_o == 1:
                    td = t + np.random.exponential(1.0 / mu)
                n_t += 1
                TravelQueueInfo.append(Na_R)
                TravelTimes.append(t + np.random.exponential(1.0 / beta))
                # update next travel event
                next_travel = min(TravelTimes) if TravelTimes else np.inf
            else:
                Utility_remote.append(0)
        # 1: local arrival
        elif idx == 1:
            ta_L = t + np.random.exponential(1.0 / lmd_L) if lmd_L > 0 else np.inf
            if n_o < nL:
                Na_L += 1
                ArrTime_L.append(t)
                n_o += 1
                OrderQueueInfo.append(-Na_L)  # negative id for local
                if n_o == 1:
                    td = t + np.random.exponential(1.0 / mu)
            else:
                Utility_local.append(0)
        # 2: service completion
        elif idx == 2:
            Nd += 1
            if t > T_w:
                cumTH += 1
            if len(OrderQueueInfo) == 0:
                # defensive check (shouldn't usually happen if td set correctly)
                td = np.inf
                continue
            n_o -= 1
            cus_id = OrderQueueInfo.pop(0)

            if cus_id > 0:  # remote completed service (not cancelled earlier)
                if cus_id not in TravelQueueInfo:
                    if t > T_w:
                        sojourn = t - ArrTime_R[cus_id - 1]
                        Utility_remote.append(V - c * sojourn)
            else:  # local
                if t > T_w:
                    sojourn = t - ArrTime_L[-cus_id - 1]
                    Utility_local.append(V - c * sojourn)

            if n_o > 0:
                td = t + np.random.exponential(1.0 / mu)
            else:
                td = np.inf

        # 3: travel completion
        else:
            # find earliest travel
            ind_t = int(np.argmin(TravelTimes))
            cus_id = TravelQueueInfo[ind_t]
            TravelTimes.pop(ind_t)
            TravelQueueInfo.pop(ind_t)

            # update next travel completion
            next_travel = min(TravelTimes) if TravelTimes else np.inf

            n_t -= 1
            if cus_id not in OrderQueueInfo:
                # customer has already been removed from order queue -> completed earlier
                if t > T_w:
                    sojourn = t - ArrTime_R[cus_id - 1]
                    Utility_remote.append(V - c * sojourn)
            else:
                QuePos = OrderQueueInfo.index(cus_id)
                # if position < nC, customer may rejoin (with prob p_R)
                if QuePos + 1 <= nC and np.random.rand() < p_R:
                    # stays in queue (nothing to do)
                    pass
                else:
                    # cancel
                    Nc += 1
                    if t > T_w:
                        sojourn = t - ArrTime_R[cus_id - 1]
                        Utility_remote.append(-c * sojourn)
                    # remove from order queue
                    OrderQueueInfo.pop(QuePos)
                    n_o -= 1
                    if n_o > 0:
                        td = t + np.random.exponential(1.0 / mu)
                    else:
                        td = np.inf

    avgUtility_remote = np.mean(Utility_remote) if len(Utility_remote) > 0 else 0.0
    avgUtility_local = np.mean(Utility_local) if len(Utility_local) > 0 else 0.0
    throughput = cumTH / (T - T_w) if (T - T_w) > 0 else 0.0
    return throughput, avgUtility_remote, avgUtility_local, mu

# --- 并行任务的Worker函数 ---
def run_simulation_for_combination(args):
    """
    为单个参数组合运行蒙特卡洛模拟，并返回所有平均指标。
    """
    current_lam, q, nR, nC, p_R, p_L, nL, V, gamma, beta, num_MC_search = args

    throughputs = np.zeros(num_MC_search)
    utils_remote = np.zeros(num_MC_search)
    utils_local = np.zeros(num_MC_search)
    service_rate_mu = 1.0

    for trial in range(num_MC_search):
        th, ur, ul, service_rate_mu = AdmissionControl(current_lam, q, p_R, p_L, nR, nC, nL, V, gamma, beta)
        throughputs[trial] = th
        utils_remote[trial] = ur
        utils_local[trial] = ul

    avg_throughput_raw = np.mean(throughputs)
    avg_throughput_capped = min(avg_throughput_raw, service_rate_mu)

    avg_utility_remote = np.mean(utils_remote)
    avg_utility_local = np.mean(utils_local)

    params_dict = {'q': q, 'p_R': p_R, 'p_L': p_L, 'nR': nR, 'nC': nC, 'nL': nL}

    return avg_throughput_raw, avg_throughput_capped, avg_utility_remote, avg_utility_local, params_dict


# --- 主程序 ---
if __name__ == "__main__":
    # --- 1. 定义4个环境参数的完整空间 ---
    # Vv = np.arange(1.5, 5.6, 1.0)        # V 的取值 [1.5, 2.5, 3.5, 4.5, 5.5]
    FIXED_V = 1.5  # <--- 在这里设定你想要的固定V值
    V_INDEX = 0  # <--- V=1.5 是第一个V值，所以它的行索引是0
    gammav = np.arange(0.1, 0.91, 0.2)   # gamma 的取值 [0.1, 0.3, 0.5, 0.7, 0.9]
    betav = np.arange(0.5, 2.51, 0.5)    # beta 的取值 [0.5, 1.0, 1.5, 2.0, 2.5]
    lmdv = np.arange(0.5, 2.51, 0.5)     # lambda 的取值 [0.5, 1.0, 1.5, 2.0, 2.5]

    # --- 2. 读取3个基准解文件中对应V_INDEX的那一行 ---
    # .iloc[V_INDEX] 会精确选取第 V_INDEX 行
    q_base_row = pd.read_csv("QE_OARC_opt2.csv", header=None).iloc[V_INDEX]
    nR_base_row = pd.read_csv("NR_OARC_opt2.csv", header=None).iloc[V_INDEX]
    nC_base_row = pd.read_csv("NC_OARC_opt2.csv", header=None).iloc[V_INDEX]

    # 将读取的行(Series)转换为包含125个值的一维Numpy数组
    q_base_flat = q_base_row.values
    nR_base_flat = nR_base_row.values
    nC_base_flat = nC_base_row.values

    # --- 3. 构建只针对固定V值的125个场景的列表 ---
    scenarios = []
    i = 0  # 计数器

    # 外层循环不再需要 V，因为它已经被固定了
    for lam in lmdv:
        for beta in betav:
            for gamma in gammav:
                scenario = {
                    "V": FIXED_V,  # <--- 使用我们定义的固定V值
                    "gamma": gamma,
                    "beta": beta,
                    "lambda": lam,
                    "q_base": q_base_flat[i],
                    "nR_base": nR_base_flat[i],
                    "nC_base": nC_base_flat[i]
                }
                scenarios.append(scenario)
                i += 1
    all_scenarios_results = []
    # addition
    p_R_base = 1.0
    p_L_base = 1.0
    nL_base = 4

    num_MC_search = 50

    print("启动针对所有125个场景的并行参数搜索...")
    start_total = time.time()

    # --- 4. 主循环：遍历我们刚刚创建的 scenarios 列表 ---
    for idx, scenario in enumerate(scenarios):
        # 从当前场景提取所有环境参数和基准解
        current_V = scenario['V']
        current_gamma = scenario['gamma']
        current_beta = scenario['beta']
        current_lam = scenario['lambda']

        q_base = scenario['q_base']
        nR_base = int(scenario['nR_base'])
        nC_base = int(scenario['nC_base'])

        print(
            f"\n--- 正在搜索场景 {idx + 1}/125: V={current_V}, gamma={current_gamma}, beta={current_beta}, lam={current_lam:.2f} ---")
        # 定义搜索空间 (逻辑不变)
        p_R_space = np.round(np.arange(0.85, 1.01, 0.05), 2)
        p_L_space = np.round(np.arange(0.85, 1.01, 0.05), 2)
        q_space = np.round(np.arange(max(0, q_base - 0.1), q_base + 0.01, 0.05), 2)
        nC_space = range(max(1, nC_base - 1), nC_base + 3)
        nR_space = range(max(1, nR_base - 2), nR_base + 3)

        # vvv --- 在这里进行修改 --- vvv
        # 根据 current_gamma 的值动态设定 nL_space 的范围
        # np.isclose 用于安全地比较浮点数
        if np.isclose(current_gamma, 0.7):
            nL_space = range(3, 7)  # gamma=0.7, nL 搜索范围是 [3, 4, 5, 6]
        elif np.isclose(current_gamma, 0.9):
            nL_space = range(3, 8)  # gamma=0.9, nL 搜索范围是 [3, 4, 5, 6, 7]
        else:
            nL_space = range(3, 6)  # 其他情况 (0.1, 0.3, 0.5), nL 搜索范围是 [3, 4, 5]
        # ^^^ ----------------------- ^^^

        param_list_for_scenario = []
        for p_R in p_R_space:
            for p_L in p_L_space:
                for q in q_space:
                    for nL in nL_space:
                        # ↓↓↓ 将 nR 的循环放在 nC 的外面 ↓↓↓
                        for nR in nR_space:
                            for nC in nC_space:
                                if nC > nR:
                                    break
                                    # vvv --- 新增剪枝逻辑 --- vvv
                                    # 判断当前组合是否被“基准解”全面压制 (Dominated)
                                if (p_R <= p_R_base and
                                        p_L <= p_L_base and
                                        q <= q_base and
                                        nL <= nL_base and
                                        nR <= nR_base and
                                        nC <= nC_base):

                                        # 如果当前组合与基准解完全相同，则不跳过，让它参与测试
                                    is_baseline = (p_R == p_R_base and p_L == p_L_base and q == q_base and
                                                    nL == nL_base and nR == nR_base and nC == nC_base)
                                    if not is_baseline:
                                        continue  # 跳过这个被压制的组合，进行下一次循环
                                    # ^^^ -------------------- ^^^

                                # 如果代码能运行到这里,
                                param_list_for_scenario.append((current_lam, q, nR, nC, p_R, p_L, nL, current_V, current_gamma, current_beta, num_MC_search))

        print(f"待测试的总组合数: {len(param_list_for_scenario)}")
        if not param_list_for_scenario:
            print("... 无有效组合可测试，跳过此场景。")
            continue

        num_processes = mp.cpu_count()
        print(f"正在将任务分配到 {num_processes} 个CPU核心...")
        with mp.Pool(processes=num_processes) as pool:
            results_for_scenario = pool.map(run_simulation_for_combination, param_list_for_scenario)

        feasible_solutions = [res for res in results_for_scenario if res[2] >= 0 and res[3] >= 0]
        print(f"  ... 找到 {len(feasible_solutions)} 个可行解。")

        # --- 5. 存储结果时，将当前场景的环境参数和最优解一起保存 ---
        scenario_result_dict = {
            'V': current_V, 'gamma': current_gamma, 'beta': current_beta, 'lambda': current_lam
        }

        if feasible_solutions:
            best_result_tuple = max(feasible_solutions, key=lambda item: item[0])
            scenario_result_dict['Max_Throughput_Raw'] = best_result_tuple[0]
            scenario_result_dict['Max_Throughput_Capped'] = best_result_tuple[1]
            scenario_result_dict['Utility_Remote'] = best_result_tuple[2]
            scenario_result_dict['Utility_Local'] = best_result_tuple[3]
            # 使用 .update() 合并两个字典
            scenario_result_dict.update(best_result_tuple[4])
        else:
            # vvv --- 在这里新增下面这行打印语句 --- vvv
            print(
                f"--- 警告: 在场景 V={current_V}, gamma={current_gamma}, beta={current_beta}, lam={current_lam:.2f} 时未找到可行解 ---")
            # ^^^ ------------------------------------- ^^^
            scenario_result_dict['Max_Throughput_Raw'] = 0
            scenario_result_dict['Max_Throughput_Capped'] = 0
            scenario_result_dict['Utility_Remote'] = -1
            scenario_result_dict['Utility_Local'] = -1
            scenario_result_dict.update({'q': -1, 'p_R': -1, 'p_L': -1, 'nR': -1, 'nC': -1, 'nL': -1})

        all_scenarios_results.append(scenario_result_dict)

    end_total = time.time()
    print(f"\n\n所有搜索完成! 总耗时: {end_total - start_total:.2f} 秒")

# --- 6. 将所有场景的结果一次性保存到Excel ---
    print("\n正在保存所有最优结果到Excel文件...")
    final_results_df = pd.DataFrame(all_scenarios_results)

# 定义最终列顺序
    column_order = ['V', 'gamma', 'beta', 'lambda',
                'Max_Throughput_Raw', 'Max_Throughput_Capped',
                'Utility_Remote', 'Utility_Local',
                'q', 'p_R', 'p_L', 'nR', 'nC', 'nL']
    final_results_df = final_results_df[column_order]

    final_results_df.to_excel("forloop_right_V1.xlsx", index=False, engine='openpyxl')
    print("最优结果已成功保存到 forloop_right_V1.xlsx")
