#!/usr/bin/env python3
"""
Script to convert markdown to HTML with MathJax support.
Uses markdown library for proper parsing.
"""

import markdown
from pathlib import Path

def markdown_to_html(md_content):
    """Convert markdown to HTML using markdown library, preserving math blocks."""
    import re
    
    # Protect math blocks by temporarily replacing them with HTML comments
    math_blocks = []
    math_counter = 0
    
    # Find and replace display math blocks ($$...$$)
    def replace_display_math(match):
        nonlocal math_counter
        math_blocks.append(match.group(0))
        placeholder = f"<!--MATH_DISPLAY_{math_counter}-->"
        math_counter += 1
        return placeholder
    
    # Find and replace inline math ($...$)
    def replace_inline_math(match):
        nonlocal math_counter
        math_blocks.append(match.group(0))
        placeholder = f"<!--MATH_INLINE_{math_counter}-->"
        math_counter += 1
        return placeholder
    
    # Protect display math ($$...$$) - non-greedy to avoid matching across blocks
    protected_content = re.sub(r'\$\$.*?\$\$', replace_display_math, md_content, flags=re.DOTALL)
    
    # Protect inline math ($...$) - but not inside display math
    # Match $...$ but not $$...$$
    protected_content = re.sub(r'(?<!\$)\$(?!\$)[^$\n]+?\$(?!\$)', replace_inline_math, protected_content)
    
    # Configure markdown with extensions
    md = markdown.Markdown(extensions=['tables', 'fenced_code'])
    html = md.convert(protected_content)
    
    # Restore math blocks (HTML comments are preserved by markdown)
    for i, math_block in enumerate(math_blocks):
        html = html.replace(f"<!--MATH_DISPLAY_{i}-->", math_block)
        html = html.replace(f"<!--MATH_INLINE_{i}-->", math_block)
    
    return html

def main():
    base_dir = Path(__file__).parent
    project_root = base_dir.parent.parent  # Go up from docs/Final_presentation to project root
    md_file = base_dir / 'methods_steps_and_results.md'
    html_file = base_dir / 'methods_steps_and_results.html'
    
    # Read markdown
    with open(md_file, 'r', encoding='utf-8') as f:
        md_content = f.read()
    
    # Convert to HTML body (basic conversion, MathJax will handle equations)
    html_body = markdown_to_html(md_content)
    
    # Fix image paths: convert relative paths (../../results/...) to paths from project root
    import re
    def fix_image_path(match):
        img_tag = match.group(0)
        # Extract the src path
        src_match = re.search(r'src="([^"]+)"', img_tag)
        if src_match:
            old_path = src_match.group(1)
            # If it's a relative path starting with ../../
            if old_path.startswith('../../'):
                # Remove the ../../ prefix to get path from project root
                new_path = old_path.replace('../../', '')
                # Verify the file exists
                abs_path = project_root / new_path
                if abs_path.exists():
                    # Use relative path from project root (works with base tag)
                    return img_tag.replace(f'src="{old_path}"', f'src="{new_path}"')
        return img_tag
    
    # Replace all img tags with fixed paths
    html_body = re.sub(r'<img[^>]+>', fix_image_path, html_body)
    
    # Wrap in HTML template with MathJax
    # Use base tag to set root for relative paths
    html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <base href="file://{project_root.absolute()}/">
    <title>Methods, Steps, and Results: Bacterial Genome Constraint Analysis</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
            line-height: 1.6;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            color: #333;
        }}
        h1, h2, h3, h4 {{
            color: #2c3e50;
            margin-top: 1.5em;
        }}
        h1 {{
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            border-bottom: 2px solid #ecf0f1;
            padding-bottom: 5px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 1em 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }}
        th {{
            background-color: #3498db;
            color: white;
        }}
        tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        img {{
            max-width: 100%;
            height: auto;
            display: block;
            margin: 20px auto;
            border: 1px solid #ddd;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        code {{
            background-color: #f4f4f4;
            padding: 2px 4px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
        }}
        pre {{
            background-color: #f4f4f4;
            padding: 10px;
            border-radius: 5px;
            overflow-x: auto;
        }}
        blockquote {{
            border-left: 4px solid #3498db;
            margin: 0;
            padding-left: 20px;
            color: #555;
        }}
    </style>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <script>
        window.MathJax = {{
            tex: {{
                inlineMath: [['$', '$']],
                displayMath: [['$$', '$$']],
                processEscapes: true,
                processEnvironments: true
            }},
            options: {{
                skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre']
            }}
        }};
    </script>
</head>
<body>
{html_body}
</body>
</html>"""
    
    # Write HTML
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(html_template)
    
    print(f"✓ Rendered {md_file.name} to {html_file.name}")

if __name__ == '__main__':
    main()

