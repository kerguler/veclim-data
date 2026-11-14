#!/bin/bash
set -e

echo "Clearing NGINX cache..."

docker compose exec -T veclim-data-server-nginx \
    sh -c 'find /var/cache/nginx -type f -delete'

docker compose exec -T veclim-data-server-nginx \
    nginx -s reload

echo "NGINX cache cleared."
