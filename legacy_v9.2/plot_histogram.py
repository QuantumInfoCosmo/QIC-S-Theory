import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# 設定：ファイルを探す候補地（これならどこにあっても見つけます）
POTENTIAL_PATHS = [
    'QIC_S_Result_N170.csv',          # 同じ場所
    'results/QIC_S_Result_N170.csv',  # resultsフォルダの中
    'data/QIC_S_Result_N170.csv'      # dataフォルダの中
]

# ファイルを自動検索
CSV_FILE = next((p for p in POTENTIAL_PATHS if os.path.exists(p)), None)
OUTPUT_IMG = 'Figure_2_Phase_Distribution.png'

def plot_histogram():
    # 1. データ読み込み
    if CSV_FILE is None:
        print("❌ エラー: CSVファイルが見つかりません！")
        print("以下の場所に 'QIC_S_Result_N170.csv' があるか確認してください：")
        for p in POTENTIAL_PATHS:
            print(f" - {p}")
        return

    print(f"✅ ファイルを発見しました: {CSV_FILE}")
    
    try:
        df = pd.read_csv(CSV_FILE)
        m_values = df['M']
    except Exception as e:
        print(f"エラー: ファイルを開けませんでした。\n{e}")
        return

    print(f"データ読み込み完了: N = {len(df)} galaxies")

    # 2. ヒストグラムの描画設定
    plt.figure(figsize=(10, 7))
    
    # ビンの設定 (論文に合わせて 0.05 刻み)
    bins = np.linspace(0, 2.0, 41)
    
    # Order相 (M < 0.5) と Chaos相 (M >= 0.5) で色分け
    order_data = m_values[m_values < 0.5]
    chaos_data = m_values[m_values >= 0.5]
    
    # 積み上げヒストグラムを描画
    plt.hist([order_data, chaos_data], bins=bins, stacked=True, 
             color=['blue', 'red'], alpha=0.7, 
             label=['Order Phase (Stable)', 'Chaos Phase (Unstable)'],
             edgecolor='black')

    # 3. デザイン調整 (論文仕様)
    plt.axvline(x=0.5, color='black', linestyle='--', linewidth=2, label='Phase Transition (M=0.5)')
    
    plt.xlabel('Phase Metric M', fontsize=14)
    plt.ylabel('Number of Galaxies', fontsize=14)
    plt.title(f'Distribution of Galactic Phases (N={len(df)})', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(axis='y', alpha=0.5)
    plt.xlim(0, 2.0)

    # 統計情報の表示 (テキストボックス)
    n_order = len(order_data)
    n_chaos = len(chaos_data)
    info_text = (f"Order Phase: {n_order} ({n_order/len(df):.1%})\n"
                 f"Chaos Phase: {n_chaos} ({n_chaos/len(df):.1%})")
    
    plt.text(1.2, plt.ylim()[1]*0.8, info_text, 
             fontsize=12, bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray'))

    # 4. 保存
    plt.tight_layout()
    plt.savefig(OUTPUT_IMG, dpi=300)
    print(f"🎉 図を保存しました: {OUTPUT_IMG}")

if __name__ == "__main__":
    plot_histogram()