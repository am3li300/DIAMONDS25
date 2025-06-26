from typing import Any, Dict
import dataclasses


@dataclasses.dataclass
class Description:
    requires_training: bool
    training_opts: Any  # Hummus should know about TrainingOptions
    hyper_params: Dict[str, Any]
