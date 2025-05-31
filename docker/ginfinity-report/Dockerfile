# Dockerfile.report – for the HTML‐report module
FROM python:3.9-slim

# Only need pandas and jinja2; no heavy deps
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
    && rm -rf /var/lib/apt/lists/*

RUN pip install pandas jinja2

WORKDIR /app
COPY nextflow/modules/generate_report.py /app/generate_report.py