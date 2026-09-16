FROM python:3.12-slim
WORKDIR /app
ENV PYTHONUNBUFFERED=1 MODELSCOPE_DATA_DIR=/data/modelscope
COPY requirements.txt requirements-live.txt ./
RUN pip install --no-cache-dir -r requirements-live.txt
COPY *.py live-settings.json ./
COPY dist/catalog.json ./dist/catalog.json
EXPOSE 8000
CMD ["python", "live_service.py", "serve"]
