import os
import dotenv
tmp = dotenv.load_dotenv()
if not tmp:
    print("Warning: .env not found!")
DIR_DATA = os.getenv("DIR_DATA")
VEC_HOST = os.getenv("VEC_HOST")
VEC_PORT = os.getenv("VEC_PORT")