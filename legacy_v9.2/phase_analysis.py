import numpy as np
import pandas as pd
import os
import glob

# 設定
DATA_DIR = 'data'
OUTPUT_CSV = 'QIC_S_Result_N170.csv'
MIN_POINTS = 5 

def calculate_properties(df):
    try:
        # データを数値に変換
        data = df.apply(pd.to_numeric, errors='coerce').dropna()
        r, v = data.iloc[:, 0].values, data.iloc[:, 1].values
        
        # 異常値の除外
        mask = (r > 0) & (v > 0)
        r, v = r[mask], v[mask]
        
        if len(r) < MIN_POINTS: return None
        
        # --- 1. Phase Metric (M) の計算 ---
        # ★ここが修正ポイント：対数(log)をとることで M=0.17 になります
        grad_H = (v ** 2) / r 
        log_grad_H = np.log(np.abs(grad_H) + 1e-10)
        m = np.var(log_grad_H)

        # --- 2. R と D_eff の計算 ---
        # ★ここも修正ポイント：この2つも計算してCSVに入れます
        r_max = np.max(r)
        v_at_rmax = v[np.argmax(r)] 
        d_eff = r_max * v_at_rmax

        return m, r_max, d_eff
    except: 
        return None

def main():
    # データフォルダ内の全ファイルを検索
    files = glob.glob(os.path.join(DATA_DIR, "*_rotmod.dat"))
    results = []
    
    print(f"Analyzing {len(files)} galaxies...")
    
    for f in files:
        try:
            # SPARCデータの読み込み
            df = pd.read_csv(f, sep=r'\s+', engine='python', header=None, comment='#')
            props = calculate_properties(df)
            if props is not None:
                m, r, d = props
                results.append({
                    'Galaxy': os.path.basename(f), 
                    'M': m,
                    'R': r,
                    'D_eff': d
                })
        except: 
            continue

    # CSVに保存（4列あることを確認！）
    res_df = pd.DataFrame(results)
    if not res_df.empty:
        # 順番を整える
        res_df = res_df[['Galaxy', 'M', 'R', 'D_eff']]
        res_df.to_csv(OUTPUT_CSV, index=False)
        print(f"Analysis Complete. Results saved to {OUTPUT_CSV}")
    else:
        print("No results found. Please check your data folder.")

if __name__ == "__main__":
    main()