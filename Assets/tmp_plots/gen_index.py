import os

def generate_index_html(directory):
    # Check if abc_tab_compare.html exists
    default_content = 'abc_tab_compare.html' if 'abc_tab_compare.html' in os.listdir(directory) else ''

    # Start of the HTML file
    html_content = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Report Index</title>
        <style>
            body { font-family: Arial, sans-serif; }
            .sidebar { 
                height: 100%;
                width: 200px;
                position: fixed;
                z-index: 1;
                top: 0;
                left: 0;
                background-color: #f1f1f1;
                overflow-x: hidden;
                padding-top: 20px;
            }
            .sidebar a {
                padding: 6px 8px 6px 16px;
                text-decoration: none;
                font-size: 16px;
                color: #818181;
                display: block;
            }
            .sidebar a:hover {
                color: #f1f1f1;
                background-color: #555;
            }
            .main { 
                margin-left: 220px;
                padding: 1px 16px;
                height: 1000px;
            }
            iframe {
                height: 100%;
                width: calc(100% - 200px);
                position: absolute;
                right: 0;
                top: 0;
                border: none;
            }
        </style>
    </head>
    <body>

    <div class="sidebar">
    <h2>Report Files</h2>
    """

    # Adding links to all HTML files in the directory
    for file in sorted(os.listdir(directory)):
        if file.endswith(".html") and file != 'index.html':
            display_name = file.replace('n_iter_', '') if file.startswith('n_iter_') else file
            html_content += f'<a href="{file}" target="reportFrame">{display_name}</a>\n'

    # End of the HTML file
    html_content += """
    </div>

    <div class="main">
        <iframe name="reportFrame" src="{default}"></iframe>
    </div>

    </body>
    </html>
    """.format(default=default_content)

    with open(os.path.join(directory, 'index.html'), 'w') as file:
        file.write(html_content)

# Replace 'hier_reports' with the path to your directory
# generate_index_html('hier_reports')
generate_index_html('reject_reports')
