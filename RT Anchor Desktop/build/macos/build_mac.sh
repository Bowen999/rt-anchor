#!/bin/bash
# build_mac.sh -- Build "RT Anchor.app" + "RT Anchor.dmg" (macOS arm64).
# Lives in build/macos/. The PyInstaller run goes to /private/tmp to avoid the
# iCloud-Desktop codesign/TCC issues; final artefacts land in dist/macos/.
#
# Requires the rt_anchor build env on the Mac (/opt/anaconda3) with rt-anchor
# installed and the pathlib backport UNINSTALLED (breaks PyInstaller).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"     # app root (build/macos/ -> ../..)
STAGE="$(mktemp -d /private/tmp/rtad_dmg.XXXXXX)"
DIST="/private/tmp/rtad_dist"
WORK="/private/tmp/rtad_build"

cd "$ROOT"

# 1. Bundle
/opt/anaconda3/bin/python -m PyInstaller --noconfirm --clean \
    --distpath "$DIST" --workpath "$WORK" \
    build/macos/RTAnchor.spec

# 2. Ad-hoc sign
codesign --force --deep --sign - "$DIST/RT Anchor.app"

# 3. Stage + compressed dmg
mkdir -p "$STAGE/Applications"
cp -R "$DIST/RT Anchor.app" "$STAGE/"
cp "$ROOT/example/samples.txt" "$ROOT/example/standards.txt" "$STAGE/" 2>/dev/null || true
hdiutil create -volname "RT Anchor" -srcfolder "$STAGE" -format UDZO \
    -ov "$ROOT/dist/macos/RT Anchor.dmg"

echo
echo "BUILD OK"
echo "  $DIST/RT Anchor.app"
echo "  $ROOT/dist/macos/RT Anchor.dmg"
echo "First open of the unsigned/ad-hoc app: right-click -> Open (or"
echo "xattr -cr \"/Applications/RT Anchor.app\") to clear Gatekeeper quarantine."
