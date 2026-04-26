import os
import argparse
import subprocess


def process_files(input_dir, quality_threshold):
    if not os.path.exists(input_dir):
        print(f"Error: input directory '{input_dir}' does not exist.")
        return

    # Get all files in the directory
    files = os.listdir(input_dir)

    # Process each file in turn
    for filename in files:
        file_path = os.path.join(input_dir, filename)

        # Check whether it is a file (exclude directories, etc.)
        if os.path.isfile(file_path):
            print(f"Processing file: {filename}")

            # Build the output filename
            output_filename = filename + ".sickle"
            output_path = os.path.join(input_dir, output_filename)

            # Build the sickle command
            command = [
                "sickle", "se",
                "-f", file_path,
                "-t", "sanger",
                "-o", output_path,
                "--quiet",
                "-q", str(quality_threshold)  # Use the user-specified quality threshold
            ]

            try:
                # Run the command
                subprocess.run(command, check=True)
                print(f"✓ Processing complete, output file: {output_filename}")
            except subprocess.CalledProcessError as e:
                print(f"✗ Processing failed: {e}")
            except Exception as e:
                print(f"✗ An error occurred: {e}")


def main():
    # Create the command-line argument parser
    parser = argparse.ArgumentParser(description="Batch process low-quality regions in FASTQ files")
    parser.add_argument("-i", "--input", required=True, help="Input folder containing FASTQ files")
    parser.add_argument("-q", "--quality", type=int, default=30,
                        help="Quality threshold parameter for sickle's -q option (default: 30)")

    # Parse arguments
    args = parser.parse_args()

    # Process files
    process_files(args.input, args.quality)


if __name__ == "__main__":
    main()
