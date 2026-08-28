#!/usr/bin/env bash
# publish_to_hisss.sh — mirror this project to the public repo github.com/CZ-Sync/HISSS
#
# Publishes exactly the git-TRACKED files of streamflowSignatures/, minus the
# exclusion list below, as a snapshot commit on HISSS main. Run it from anywhere;
# it locates itself. Re-run after merging work to keep the public mirror current.
#
# Exclusions (decided 2026-08-28, pre-publication audit):
#   - Claude tooling            (CLAUDE.md, claude-skill/, .claudeignore)
#   - Unpublished manuscript    (docs/MANUSCRIPT_DRAFT.md)
#   - Internal planning records (docs/plans/)
#   - Legacy Shiny app          (streamflowAndClimateVisualizationApp/ — defunct S3 backend)
#   - Large validation CSVs     (golden-outputs/*.csv > 50 MB; README + metadata CSV stay)
#   - This is a MIRROR: do not commit to HISSS directly — changes there will be
#     overwritten on the next publish.

set -euo pipefail

SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REMOTE="git@github.com:CZ-Sync/HISSS.git"
BRANCH="main"

EXCLUDE_PATTERNS=(
  '^CLAUDE\.md$'
  '^\.claudeignore$'
  '^claude-skill/'
  '^docs/MANUSCRIPT_DRAFT\.md$'
  '^docs/plans/'
  '^streamflowAndClimateVisualizationApp/'
  '^golden-outputs/streamflow_signatures_julia_apr2026\.csv$'
  '^golden-outputs/streamflow_signatures_full_10feb2026\.csv$'
)

WORK="$(mktemp -d /tmp/hisss_publish.XXXXXX)"
trap 'rm -rf "$WORK"' EXIT

echo "Cloning $REMOTE ..."
git clone --quiet --depth 1 "$REMOTE" "$WORK" 2>/dev/null || {
  # empty repo: clone still succeeds but warns; init fallback for safety
  git -C "$WORK" init --quiet 2>/dev/null || true
}
git -C "$WORK" checkout -q -B "$BRANCH" 2>/dev/null || true

# Wipe everything except .git so deletions in the source propagate
find "$WORK" -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +

# Copy exactly the tracked files, minus exclusions
cd "$SRC"
FILTER="$(printf '%s|' "${EXCLUDE_PATTERNS[@]}")"
FILTER="${FILTER%|}"
git ls-files | grep -Ev "$FILTER" | while IFS= read -r f; do
  mkdir -p "$WORK/$(dirname "$f")"
  cp -p "$f" "$WORK/$f"
done

SRC_REV="$(git rev-parse --short HEAD)"
DIRTY=""
git diff --quiet HEAD -- . 2>/dev/null || DIRTY=" (+uncommitted changes)"

cd "$WORK"
# -f: the mirrored .gitignore must not suppress files tracked in the source repo
# (e.g. streamflowSignatures.Rproj is tracked there despite the *.Rproj rule)
git add -A -f
if git diff --cached --quiet; then
  echo "HISSS is already up to date with $SRC_REV — nothing to publish."
  exit 0
fi
git -c user.name="$(git -C "$SRC" config user.name)" \
    -c user.email="$(git -C "$SRC" config user.email)" \
    commit --quiet -m "Sync from source repo @ ${SRC_REV}${DIRTY} ($(date +%Y-%m-%d))"
git push --quiet "$REMOTE" "$BRANCH"
echo "Published to https://github.com/CZ-Sync/HISSS @ source ${SRC_REV}${DIRTY}"
