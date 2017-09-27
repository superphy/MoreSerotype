import logging
import os

import definitions

logger = logging.getLogger(__name__)

def setup_logging():
    # set up logging to file - see previous section for more details
    if not os.path.exists(definitions.TEMP_DIR):
        os.makedirs(definitions.TEMP_DIR)
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m/%d %H:%M',
                        filename='temp/main.log',
                        filemode='w')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter(
        fmt='%(name)-12s: %(levelname)-8s %(message)s',
        datefmt='%H:%M:%S'
    )
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # Now, we can log to the root logger, or any other logger. First the root...
    logger.info('Logger initialized.')
