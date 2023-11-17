#!/bin/sh

gh api -H "Accept: application/vnd.github+json" -H "X-GitHub-Api-Version: 2022-11-28" /repos/kkmaute/moris/traffic/clones > /home/maute/codes/moris/share/maintenance/github_traffic/`date +%d-%m-%y`
