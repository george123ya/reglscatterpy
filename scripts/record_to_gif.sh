#!/usr/bin/env bash
# Convert a screen recording (mp4/webm/mkv) into a small, high-quality GIF for the
# README — two-pass palette so colours stay crisp at a few hundred KB.
#
#   scripts/record_to_gif.sh in.mp4 assets/demo.gif [WIDTH] [FPS] [START] [DUR]
#
# Examples:
#   scripts/record_to_gif.sh ~/Videos/cap.mp4 assets/demo.gif            # 760px, 15fps
#   scripts/record_to_gif.sh ~/Videos/cap.mp4 assets/demo.gif 760 15 3 9 # trim 3s..12s
#
# Needs ffmpeg. Record the source with whatever you already use (OBS, the omarchy
# screen-record hotkey, etc.); this just makes the optimised GIF.
set -euo pipefail

IN="${1:?usage: record_to_gif.sh in.mp4 out.gif [width] [fps] [start] [dur]}"
OUT="${2:?missing output path, e.g. assets/demo.gif}"
WIDTH="${3:-760}"
FPS="${4:-15}"
START="${5:-}"
DUR="${6:-}"

trim=()
[[ -n "$START" ]] && trim+=(-ss "$START")
[[ -n "$DUR" ]] && trim+=(-t "$DUR")

PAL="$(mktemp --suffix=.png)"
trap 'rm -f "$PAL"' EXIT

filters="fps=${FPS},scale=${WIDTH}:-1:flags=lanczos"
ffmpeg -hide_banner -loglevel error -y "${trim[@]}" -i "$IN" \
    -vf "${filters},palettegen=stats_mode=diff" "$PAL"
ffmpeg -hide_banner -loglevel error -y "${trim[@]}" -i "$IN" -i "$PAL" \
    -lavfi "${filters}[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=3" "$OUT"

bytes=$(stat -c%s "$OUT")
printf 'wrote %s  (%d KB, %dpx @ %dfps)\n' "$OUT" "$((bytes / 1024))" "$WIDTH" "$FPS"
