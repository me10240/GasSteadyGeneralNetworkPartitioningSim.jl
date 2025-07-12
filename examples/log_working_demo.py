
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG) # decides which events should be propagated
import network_partition_script as nps

class CustomLogFormatter(logging.Formatter):
    format = "%(asctime)s - %(name)s  - %(levelname)s - (%(filename)s:%(lineno)d  %(funcName)s) - %(message)s"
    FORMATS = {
        logging.DEBUG: format,
        logging.INFO: format,
        logging.WARNING: format,
        logging.ERROR: format,
        logging.CRITICAL: format
    }
    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)


class CustomStreamFormatter(logging.Formatter):

    green = "\x1b[32;20m"
    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(name)s - %(levelname)s - (%(filename)s:%(lineno)d  %(funcName)s) - %(message)s "

    FORMATS = {
        logging.DEBUG: green + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt, "%Y-%m-%d %H:%M:%S")
        return formatter.format(record)




def main():
    
    # this is how you access the logger of a module you are using and turn it off (or not)
    nps_logger = logging.getLogger(nps.__name__) # works because named it with good programming practice
    # but what if inside module, a different  name was hardwired ? then wont work unless that name is used
    # what if there is no logger ? will this fail gracefully ? 
    nps_logger.disabled = True


    dirname = "./data/Texas7k_Gas/"
    loglevel_str = "info"

    nps.run_script(dirname, loglevel="info", allow_slack_node_partitioning = False, num_max=2, round_max=20, plotting_flag=False)
    
    import os
    level = getattr(logging, loglevel_str.upper())  #levels are 10, 20, 30, 40, 50

    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)
    ch.setFormatter(CustomStreamFormatter()) 
    log.addHandler(ch)
    
    # create file handler
    logfile = dirname + "dummy.log"
    fh = logging.FileHandler(logfile, mode='w')
    fh.setLevel(logging.DEBUG) 
    fh.setFormatter(CustomLogFormatter())
    log.addHandler(fh) 

    log.info(os.getcwd())

    log.info("Using Logger version {}".format(logging.__version__))
    # dirname = "8-node/"
    

    
if __name__ == "__main__":
    main()