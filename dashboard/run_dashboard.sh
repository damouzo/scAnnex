#!/bin/bash

# scAnnex Dashboard Runner
# Builds and runs the Docker-based Shiny dashboard

set -euo pipefail

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Configuration
IMAGE_NAME="scannex-dashboard"
TAG="latest"
CONTAINER_NAME="scannex-dashboard"
PORT="3838"

# Directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
TEST_DATA_DIR="$PROJECT_ROOT/test_data/analytical_core_results"

echo "======================================================================="
echo "scAnnex Dashboard - Docker Runner"
echo "======================================================================="
echo ""

# ==============================================================================
# Parse command line arguments
# ==============================================================================

ACTION="${1:-run}"

case "$ACTION" in
  
  build)
    echo "Building Docker image..."
    docker build -t "$IMAGE_NAME:$TAG" "$SCRIPT_DIR"
    echo -e "\n${GREEN}✓ Docker image built successfully${NC}"
    echo "  Image: $IMAGE_NAME:$TAG"
    ;;
    
  run)
    echo "Starting scAnnex Dashboard..."
    echo ""
    
    # Check if test data exists
    if [ ! -d "$TEST_DATA_DIR" ]; then
      echo -e "${YELLOW}⚠ Test data not found: $TEST_DATA_DIR${NC}"
      echo "  Please run test_data/test_analytical_core.sh first to generate test data"
      exit 1
    fi
    
    # Check if container already running
    if docker ps -a --format '{{.Names}}' | grep -q "^${CONTAINER_NAME}$"; then
      echo "Stopping existing container..."
      docker stop "$CONTAINER_NAME" || true
      docker rm "$CONTAINER_NAME" || true
    fi
    
    # Check if image exists, build if needed
    if ! docker images --format '{{.Repository}}:{{.Tag}}' | grep -q "^${IMAGE_NAME}:${TAG}$"; then
      echo "Image not found, building..."
      docker build -t "$IMAGE_NAME:$TAG" "$SCRIPT_DIR"
    fi
    
    echo "Launching dashboard container..."
    docker run -d \
      -p "$PORT:3838" \
      --name "$CONTAINER_NAME" \
      -v "$TEST_DATA_DIR:/srv/shiny-server/data" \
      "$IMAGE_NAME:$TAG"
    
    echo -e "\n${GREEN}✓ Dashboard started successfully${NC}"
    echo ""
    echo "Access the dashboard at:"
    echo "  http://localhost:$PORT"
    echo ""
    echo "To view logs:"
    echo "  docker logs -f $CONTAINER_NAME"
    echo ""
    echo "To stop:"
    echo "  $0 stop"
    ;;
    
  stop)
    echo "Stopping dashboard..."
    docker stop "$CONTAINER_NAME" || true
    docker rm "$CONTAINER_NAME" || true
    echo -e "${GREEN}✓ Dashboard stopped${NC}"
    ;;
    
  restart)
    "$0" stop
    sleep 2
    "$0" run
    ;;
    
  logs)
    docker logs -f "$CONTAINER_NAME"
    ;;
    
  shell)
    echo "Opening shell in container..."
    docker exec -it "$CONTAINER_NAME" /bin/bash
    ;;
    
  *)
    echo "Usage: $0 {build|run|stop|restart|logs|shell}"
    echo ""
    echo "Commands:"
    echo "  build    - Build the Docker image"
    echo "  run      - Build (if needed) and run the dashboard"
    echo "  stop     - Stop the dashboard container"
    echo "  restart  - Restart the dashboard"
    echo "  logs     - Show dashboard logs (follow mode)"
    echo "  shell    - Open a shell in the running container"
    exit 1
    ;;
esac
