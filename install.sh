#!/bin/sh
# Foci Analyzer: installation on macOS and Linux. In a terminal, in the folder of the program:
#     sh install.sh
# Run it again to update or to add networks. Everything goes into one folder (default ~/FociAnalyzer, or
# $FOCI_ANALYZER_INSTALL_DIR); removing it (and the app or menu entry) removes the program.
# No Python is needed beforehand: uv (https://docs.astral.sh/uv/) fetches it. Files fetched like this are not
# marked as downloaded by macOS, so the program starts without a warning, also unsigned.
# Options are passed on to installer/install.py, for example:  sh install.sh --networks stardist,spotiflow
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
TARGET=${FOCI_ANALYZER_INSTALL_DIR:-$HOME/FociAnalyzer}
echo
echo " Foci Analyzer will be installed in $TARGET"
echo
mkdir -p "$TARGET/bin"
UV="$TARGET/bin/uv"
if [ ! -x "$UV" ]; then
    echo " Fetching uv (a small tool that installs Python and the packages)..."
    if command -v curl >/dev/null 2>&1; then
        curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR="$TARGET/bin" UV_NO_MODIFY_PATH=1 sh
    else
        wget -qO- https://astral.sh/uv/install.sh | env UV_INSTALL_DIR="$TARGET/bin" UV_NO_MODIFY_PATH=1 sh
    fi
fi
export UV_PYTHON_INSTALL_DIR="$TARGET/python" UV_CACHE_DIR="$TARGET/cache"
"$UV" run --no-project --python 3.12 "$HERE/installer/install.py" --target "$TARGET" --uv "$UV" --shortcut "$@"
