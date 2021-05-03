import logging

def get_logger(log_file_path):
# Create a custom logger
    logger = logging.getLogger(__name__)

# Create handlers
    stream_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(log_file_path)
    stream_handler.setLevel(logging.INFO)
    file_handler.setLevel(logging.ERROR)

    # Create formatters and add it to handlers
    stream_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    file_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(stream_format)
    file_handler.setFormatter(file_format)

    # Add handlers to the logger
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    return logger
