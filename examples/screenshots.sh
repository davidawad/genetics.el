#!/bin/sh
# Regenerate docs/screenshots/*.png from the synthetic test fixtures.
#
#   examples/screenshots.sh
#
# Needs a graphical Emacs (EMACS, default `emacs'), xwd and ImageMagick
# (`magick' or `convert').  When DISPLAY is unset it starts Xvfb (XVFB,
# default `Xvfb') on :99 for the run.  Debian: apt install xvfb x11-apps
# imagemagick emacs-gtk fonts-jetbrains-mono.  Only synthetic data is
# opened (test/fixtures), with `emacs -Q'.
set -eu
root=$(cd "$(dirname "$0")/.." && pwd)
EMACS=${EMACS:-emacs}
XVFB=${XVFB:-Xvfb}

# Keep the sample report in sync with the code before photographing it.
"$EMACS" -Q --batch -l "$root/examples/regenerate-report.el"

xvfb_pid=
if [ -z "${DISPLAY:-}" ]; then
  "$XVFB" :99 -screen 0 1600x1000x24 -nolisten tcp >/dev/null 2>&1 &
  xvfb_pid=$!
  trap 'kill $xvfb_pid 2>/dev/null' EXIT
  export DISPLAY=:99
  sleep 2
fi

cd "$root"
"$EMACS" -Q -l "$root/examples/screenshots.el"
ls -l "$root"/docs/screenshots/*.png
