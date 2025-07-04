from pathlib import Path
import os
from dotenv import load_dotenv

load_dotenv()

NCIT_API_KEY = os.getenv("NCIT_API_KEY")
NCBI_API_KEY = os.getenv("NCBI_API_KEY")
EMAIL = os.getenv("EMAIL_ADDRESS")
UMLS_API_KEY = os.getenv("UMLS_API_KEY")
