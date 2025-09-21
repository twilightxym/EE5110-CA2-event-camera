#!/usr/bin/env bash
set -euo pipefail

# 处理当前目录（及其子目录）下的所有视频
INPUT_DIR="."
FPS_CFR=""   # 可留空；如需强制容器帧率可写：-r 30

# 检查 ffmpeg
if ! command -v ffmpeg >/dev/null 2>&1; then
  echo "Error: ffmpeg 未安装或不在 PATH。请先安装 ffmpeg。"
  exit 1
fi

# 递归查找视频文件
find "$INPUT_DIR" -type f \( -iname "*.mp4" -o -iname "*.mov" -o -iname "*.m4v" \) | \
while IFS= read -r inpath; do
  indir="$(dirname "$inpath")"     # 原视频所在目录（比如 ./encoded_test/Type1）
  base="$(basename "$inpath")"     # 文件名含后缀
  name="${base%.*}"                # 去后缀的文件名

  outdir="${indir}/2K"             # 在同级建 2K 文件夹
  outfile="${outdir}/${name}_2k.mp4"

  # 已存在就跳过（如需覆盖，删掉这个 if 块）
  if [[ -f "$outfile" ]]; then
    echo "Skip (exists): $outfile"
    continue
  fi

  mkdir -p "$outdir"

  echo "Processing:"
  echo "  IN : $inpath"
  echo "  OUT: $outfile"

  # 由于原视频是 DCI 4K(4096×2160)，严格 1/2 降采样到 2048×1080（整数倍，不会破坏比例）
  # 若要容错其它比例，可改成: -vf "scale=2048:-2"
  ffmpeg -y -loglevel error -i "$inpath" \
    -vf "scale=2048:1080" \
    -c:v libx264 -pix_fmt yuv420p -vsync cfr ${FPS_CFR} \
    -an \
    "$outfile"

  echo "Done."
done

echo "All finished. Each subfolder now has its own '2K' with outputs."