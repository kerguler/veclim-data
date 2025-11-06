import os
import dotenv
tmp = dotenv.load_dotenv()
if not tmp:
    print("Warning: .env not found!")
DIR_DATA  = os.getenv("DIR_DATA")
DIR_CACHE = os.getenv("DIR_CACHE")
DIR_TILE  = os.getenv("DIR_TILE")
VEC_HOST  = os.getenv("VEC_HOST")