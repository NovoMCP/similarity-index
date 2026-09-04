# NovoMCP similarity index — FPSim2 exact-Tanimoto search over the open corpus.
FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY build_index.py server.py entrypoint.sh ./
RUN chmod +x entrypoint.sh

# Serve mode defaults; build mode overrides via `docker run ... build ...`.
ENV INDEX_DIR=/data/index \
    INDEX_MODE=in-memory \
    PORT=8080
EXPOSE 8080

ENTRYPOINT ["./entrypoint.sh"]
CMD ["serve", "--index", "/data/index"]
