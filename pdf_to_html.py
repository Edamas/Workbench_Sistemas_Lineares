import os
from pathlib import Path
import fitz # PyMuPDF
import glob

# --- Configuration ---
PDF_DIRS = {
    "Notas": Path("D:/UNIVESP/Qualificações profissionais/USP/Álgebra Linear - IME/notas"),
    "Listas": Path("D:/UNIVESP/Qualificações profissionais/USP/Álgebra Linear - IME/listas"),
    "Provas": Path("D:/UNIVESP/Qualificações profissionais/USP/Álgebra Linear - IME/provas"),
}
HTML_OUTPUT_DIR = Path("html_output")
HTML_OUTPUT_DIR.mkdir(exist_ok=True)

# --- Functions ---

def get_all_pdf_files():
    """Scans all configured directories and returns a flat list of PDF files."""
    all_files = []
    for category, path in PDF_DIRS.items():
        files = sorted(path.glob("*.pdf"))
        all_files.extend(files)
    return all_files

def extract_html_from_pdf(pdf_path: Path, output_dir: Path):
    """
    Extracts content from a PDF file as HTML and saves it to a file.
    """
    html_content = ""
    try:
        document = fitz.open(pdf_path)
        for page_num in range(document.page_count):
            page = document.load_page(page_num)
            html_content += page.get_text("html")
        document.close()
        
        output_path = output_dir / f"{pdf_path.stem}.html"
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html_content)
        print(f"Successfully extracted HTML from '{pdf_path.name}' to '{output_path.name}'")

    except Exception as e:
        print(f"Error extracting HTML from {pdf_path}: {e}")

# --- Main Script ---

if __name__ == "__main__":
    all_pdfs = get_all_pdf_files()
    
    if not all_pdfs:
        print("No PDF files found in the configured directories.")
    else:
        print(f"Found {len(all_pdfs)} PDF files to process.")
        for pdf_path in all_pdfs:
            extract_html_from_pdf(pdf_path, HTML_OUTPUT_DIR)
        print("\nHTML extraction complete.")
