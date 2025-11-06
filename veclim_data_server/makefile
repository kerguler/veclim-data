PROJECT_NAME := $(shell basename "$$(pwd)")

.PHONY: deploy down up rebuild clean-cache logs

# Stop stack and remove nginx cache volume, then rebuild & start
deploy: down clean-cache rebuild up

# Stop containers
down:
	@echo "Stopping containers..."
	docker compose down

# Remove only the nginx cache volume
clean-cache:
	@echo "Removing NGINX cache volume..."
	-volume_name=$$(docker volume ls -q | grep -E "$(PROJECT_NAME)_nginx_cache"); \
	if [ -n "$$volume_name" ]; then \
		echo "Deleting cache volume: $$volume_name"; \
		docker volume rm $$volume_name; \
	else \
		echo "No nginx cache volume found (ok)"; \
	fi

# Rebuild images
rebuild:
	@echo "Rebuilding images..."
	docker compose build

# Start stack
up:
	@echo "Starting containers..."
	docker compose up -d

# View logs
logs:
	docker compose logs -f
