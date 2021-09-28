# Mapping Teichlab computational cloud

A set of scripts for mapping 10X data:
  - Cellranger supported via `scripts/pyranger.py`
  - Starsolo supported via `runscripts/starsolo.sh`
  - Kraken metagenomics supported via `runscripts/kraken.sh`
  - Visium supported via `runscripts/visium.sh`

The shell based scripts assume this repo present in `/mnt/mapcloud`, pyranger infers the location itself and is less picky.

On the off chance you read this and care about something, contact kp9 on Mattermost for details.