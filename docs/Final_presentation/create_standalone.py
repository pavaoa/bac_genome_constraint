#!/usr/bin/env python3
"""
Create a standalone HTML package with all images included.
This allows the HTML to be viewed without access to the full repository.
"""

import re
import shutil
from pathlib import Path
from zipfile import ZipFile

def create_standalone():
    """Create a standalone HTML package with all images."""
    base_dir = Path(__file__).parent
    project_root = base_dir.parent.parent
    html_file = base_dir / 'methods_steps_and_results.html'
    standalone_dir = base_dir / 'standalone'
    images_dir = standalone_dir / 'images'
    
    # Create directories
    standalone_dir.mkdir(exist_ok=True)
    images_dir.mkdir(exist_ok=True)
    
    print(f"Creating standalone package in {standalone_dir}...")
    
    # Read the HTML file
    with open(html_file, 'r', encoding='utf-8') as f:
        html_content = f.read()
    
    # Find all image references
    image_pattern = r'src="([^"]+\.(?:png|pdf|jpg|jpeg|gif|svg))"'
    image_matches = re.findall(image_pattern, html_content, re.IGNORECASE)
    
    print(f"Found {len(image_matches)} image references")
    
    # Track copied images to avoid duplicates
    copied_images = {}
    image_counter = 1
    
    # Copy images and update paths
    def replace_image_path(match):
        nonlocal image_counter
        full_match = match.group(0)
        img_path = match.group(1)
        
        # Skip if already processed
        if img_path in copied_images:
            new_path = copied_images[img_path]
            return full_match.replace(img_path, new_path)
        
        # Determine source path
        if img_path.startswith('results/'):
            # Path relative to project root
            source_path = project_root / img_path
        elif img_path.startswith('file://'):
            # Absolute file:// path
            source_path = Path(img_path.replace('file://', ''))
        elif img_path.startswith('../'):
            # Relative path from HTML location
            source_path = base_dir / img_path
        else:
            # Assume relative to project root
            source_path = project_root / img_path
        
        if not source_path.exists():
            print(f"  ⚠ Warning: Image not found: {source_path}")
            return full_match
        
        # Create a clean filename
        # Extract the original filename
        original_name = source_path.name
        # Create a unique name if needed (handle duplicates)
        if original_name in [v for v in copied_images.values()]:
            name_parts = original_name.rsplit('.', 1)
            if len(name_parts) == 2:
                new_name = f"{name_parts[0]}_{image_counter}.{name_parts[1]}"
            else:
                new_name = f"{original_name}_{image_counter}"
        else:
            new_name = original_name
        
        # Copy image to images directory
        dest_path = images_dir / new_name
        try:
            shutil.copy2(source_path, dest_path)
            print(f"  ✓ Copied: {source_path.name} -> images/{new_name}")
        except Exception as e:
            print(f"  ✗ Error copying {source_path}: {e}")
            return full_match
        
        # Store mapping and update path
        new_relative_path = f"images/{new_name}"
        copied_images[img_path] = new_relative_path
        image_counter += 1
        
        return full_match.replace(img_path, new_relative_path)
    
    # Replace all image paths
    updated_html = re.sub(image_pattern, replace_image_path, html_content, flags=re.IGNORECASE)
    
    # Remove the base tag since we're using relative paths now
    updated_html = re.sub(r'<base[^>]*>', '', updated_html)
    
    # Fix math equation rendering: Ensure math blocks are properly formatted for MathJax
    # MathJax needs math blocks to be in block-level elements, not wrapped in <p> tags
    # Pattern 1: Remove <p> tags that directly wrap math blocks
    updated_html = re.sub(r'<p>\s*(\$\$.*?\$\$)\s*</p>', r'\n\1\n', updated_html, flags=re.DOTALL)
    
    # Pattern 2: Wrap standalone math blocks in divs for better MathJax processing
    # Find math blocks that are on their own lines (not in tags)
    def wrap_math_block(match):
        math_content = match.group(1)
        # Clean up whitespace
        math_content = math_content.strip()
        return f'\n<div class="math-display">\n{math_content}\n</div>\n'
    
    # Match math blocks that are standalone (preceded by > or newline, followed by < or newline)
    updated_html = re.sub(r'(?<=[>\n])\s*(\$\$[^$]+\$\$)\s*(?=[<\n])', wrap_math_block, updated_html, flags=re.DOTALL)
    
    # Clean up any excessive whitespace
    updated_html = re.sub(r'\n{3,}', '\n\n', updated_html)
    
    # Add CSS for math-display divs to ensure proper spacing
    math_display_css = """
        .math-display {
            margin: 1.5em 0;
            text-align: center;
            overflow-x: auto;
        }
    """
    # Insert the CSS before </style>
    if '</style>' in updated_html:
        updated_html = updated_html.replace('</style>', f'{math_display_css}\n    </style>')
    
    # Write updated HTML
    standalone_html = standalone_dir / 'methods_steps_and_results.html'
    with open(standalone_html, 'w', encoding='utf-8') as f:
        f.write(updated_html)
    
    print(f"✓ Created standalone HTML: {standalone_html}")
    print(f"✓ Copied {len(copied_images)} unique images to {images_dir}")
    
    # Create zip file
    zip_file = base_dir / 'methods_steps_and_results_standalone.zip'
    print(f"\nCreating zip file: {zip_file}")
    
    with ZipFile(zip_file, 'w') as zipf:
        # Add HTML file
        zipf.write(standalone_html, standalone_html.name)
        print(f"  ✓ Added {standalone_html.name}")
        
        # Add all images
        for img_file in images_dir.glob('*'):
            if img_file.is_file():
                zipf.write(img_file, f"images/{img_file.name}")
                print(f"  ✓ Added images/{img_file.name}")
    
    print(f"\n✓ Standalone package created successfully!")
    print(f"  HTML: {standalone_html}")
    print(f"  Images: {images_dir} ({len(list(images_dir.glob('*')))} files)")
    print(f"  Zip: {zip_file}")
    print(f"\nTo use: Extract the zip file and open methods_steps_and_results.html in a browser.")

if __name__ == '__main__':
    create_standalone()

