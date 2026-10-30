# 成功したGaussian固定断片SPの保存

4系（RAFA14・MoriPF・RAFA15・RAFA16）×3点×2断片、計24件。各フォルダに実際の入力 `fragment.gjf` と成功した出力全文 `fragment.out.gz` を収録する。Gaussianソフト本体、chk/rwf、失敗・中断試行は含めない。

`index.json` に入力・gzip・展開後出力のSHA-256、元出力サイズを記録。各系の `manifest.json` は座標・原子対応、既存全系エネルギーと分離反応物参照、元作業フォルダの来歴を保存する。Windows形式の元パスは来歴であり、別環境でそのパスの存在を要求しない。元全系・参照ログは既存 `data/completed_20261017/` の追跡済みデータを参照する。

リポジトリルートで、Python環境に `irc_eda/requirements.txt` の依存関係を導入した後：

```text
python -X utf8 irc_eda/scripts/archive_smd_pilots.py verify
python -X utf8 irc_eda/scripts/summarize_smd_all_pilots.py
```

1行目は24件を展開・hash照合・電子条件と座標の検証、2行目はさらに各系CSV/図/報告と4系比較を再生成する。両方ともGaussianを起動しない。Windows固有パスや元outputsフォルダに依存しない。入力は `.gitattributes` で改行変換を無効化し、圧縮出力は元バイト列を保存する。
