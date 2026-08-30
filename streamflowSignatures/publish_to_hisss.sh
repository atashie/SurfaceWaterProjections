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

# Distinguish "remote is empty" (init locally) from real failures (abort loudly):
# ls-remote --exit-code returns 2 for a reachable repo with no refs, 128 for
# auth/DNS/not-found — never treat those as an empty repo.
echo "Checking $REMOTE ..."
if git ls-remote --exit-code "$REMOTE" >/dev/null; then
  git clone --quiet --depth 1 "$REMOTE" "$WORK"
else
  rc=$?
  if [ "$rc" -eq 2 ]; then
    echo "Remote is empty — initializing first publish."
    git init --quiet "$WORK"
  else
    echo "ERROR: cannot reach $REMOTE (git ls-remote exit $rc) — aborting." >&2
    exit 1
  fi
fi
git -C "$WORK" checkout -q -B "$BRANCH"

# Wipe everything except .git so deletions in the source propagate
find "$WORK" -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +

# Copy exactly the tracked files, minus exclusions
cd "$SRC"
FILTER="$(printf '%s|' "${EXCLUDE_PATTERNS[@]}")"
FILTER="${FILTER%|}"
# Collect the list first so an empty selection fails deliberately (under
# pipefail a no-match grep would otherwise abort the script with no message).
FILES="$(git ls-files | grep -Ev "$FILTER" || true)"
if [ -z "$FILES" ]; then
  echo "ERROR: publish set is empty after exclusions — refusing to publish." >&2
  exit 1
fi
printf '%s\n' "$FILES" | while IFS= read -r f; do
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
