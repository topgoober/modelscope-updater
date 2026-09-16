# Modelscope updater

Python ingestion service for NOAA and ECMWF forecasts. Supplies read-only forecast frames to the existing Modelscope map. Uses a Dockerfile; start command is `python live_service.py serve`.

## Hosting

Attach a persistent volume at `/data`. Set `MODELSCOPE_DATA_DIR=/data/modelscope`, `MODELSCOPE_EPHEMERAL_HISTORY=1`, `OPENBLAS_NUM_THREADS=1`, and `OMP_NUM_THREADS=1`. The service listens on the platform-provided `PORT` (default 8000). Health check: `/healthz`. Use one replica and disable sleeping.

The feed serves public meteorological data. CORS is restricted by `MODELSCOPE_SITE_ORIGIN`, defaulting to `https://modelscope-weather.wilcob9139.chatgpt.site`. There are no upload or administrative endpoints.

## Forecast coverage

Nine model families, selected EPS/SREF members, 3-hour valid times through the next 48 hours. Available fields and profiles depend on upstream data. History-derived QPF/UH require complete source history. Missing data is never invented. RRFS uses the configured parallel feed.

New inventories are checked on a 10-minute idle/retry interval. Processing adds latency; this is not a guarantee of delivery within 10 minutes. Browser refresh is separately configured on the website. A healthy service does not imply that every model has arrived.

`MODELSCOPE_SETTINGS` accepts JSON overrides for `live-settings.json`. Published frames default to a 320 MB storage budget, with distant valid times removed first if required. Per-job raw caches are deleted after processing. The code is prepared for a small worker, but actual memory use, ingestion speed, and source availability must be verified after deployment.

The catalog bundled under `dist/` supplies parameter definitions; example frame files are not included. The service begins with an empty live catalog and adds real forecasts as downloads succeed. The website must be pointed to the deployed feed only after ingestion is verified.

## Tests

Install `requirements-live.txt`, then run `python -m unittest tests.test_live`. Tests cover orchestration and read-only serving; they do not substitute for live ingestion verification.
